#=
    Designer based on the equivalent annual cost (EAC) with multiple scenarios
=#

mutable struct MILPOptions
  solver::Module
  reducer::AbstractScenariosReducer
  objective_risk::AbstractRiskMeasure
  share_risk::AbstractRiskMeasure
  reopt::Bool
  read_reduction::Union{String, Nothing}
  write_reduction::Union{String, Nothing}

  MILPOptions(; solver = CPLEX,
                reducer = FeatureBasedReducer(),
                objective_risk = Expectation(),
                share_risk = Expectation(),
                reopt = false,
                read_reduction = nothing,
                write_reduction = nothing) =
                new(solver, reducer, objective_risk, share_risk, reopt, read_reduction, write_reduction)
end

mutable struct MILP <: AbstractDesigner
    options::MILPOptions
    decisions::NamedTuple
    model::JuMP.Model
    history::AbstractScenarios

    MILP(; options = MILPOptions()) = new(options)
end

### Models
function build_model(mg::Microgrid, designer::MILP, ω::Scenarios, probabilities::Vector{Float64})
    # Initialize
    m = Model(designer.options.solver.Optimizer)
    # Add decision variables
    add_decisions!(m, mg, ω)
    # Add technical constraints
    add_technical_constraints!(m, mg, ω)
    # Add power balance constraints
    add_power_balances!(m, mg, ω)
    # Renewable share constraint
    add_renewable_share!(m, mg, ω, probabilities, designer.options.share_risk)
    # Objective
    add_objective!(m, mg, ω, probabilities, designer.options.objective_risk)
    return m
end

# Decisions
function add_decisions!(m::Model, mg::Microgrid, ω::Scenarios)
    nh, ns = size(ω.demands[1].power,1), size(ω.demands[1].power,3)
    # Generation
    if !isempty(mg.generations)
        na = length(mg.generations)
        # Investment variables
        @variable(m, r_g[1:na])
    end
    # Storages
    if !isempty(mg.storages)
        na = length(mg.storages)
        # Investment variables
        @variable(m, r_sto[1:na])
        # Operation variables
        @variables(m, begin
        p_ch[1:nh, 1:ns, 1:na]   >= 0.
        p_dch[1:nh, 1:ns, 1:na]  >= 0.
        soc[1:nh+1, 1:ns, 1:na]
        end)
    end
    # Converters
    if !isempty(mg.converters)
        na = length(mg.converters)
        # Investment variables
        @variable(m, r_c[1:na])
        @variable(m, p_c[1:nh, 1:ns, 1:na] >= 0.)
    end
    # Grids
    if !isempty(mg.grids)
        na = length(mg.grids)
        # Operation variables
        @variables(m, begin
        p_in[1:nh, 1:ns, 1:na]   >= 0.
        p_out[1:nh, 1:ns, 1:na]  >= 0.
        end)
    end
end
# Technical constraints
function add_technical_constraints!(m::Model, mg::Microgrid, ω::Scenarios)
    nh, ns = size(ω.demands[1].power,1), size(ω.demands[1].power,3)
    # Generations
    if !isempty(mg.generations)
        na = length(mg.generations)
        # Constraint
        @constraints(m, begin
        # Investment bounds
        [a in 1:na], m[:r_g][a] >= mg.generations[a].bounds.lb
        [a in 1:na], m[:r_g][a] <= mg.generations[a].bounds.ub
        end)
    end
    # Storages
    if !isempty(mg.storages)
        na = length(mg.storages)
        # Constraint
        @constraints(m, begin
        # Investment bounds
        [a in 1:na], m[:r_sto][a] >= mg.storages[a].bounds.lb
        [a in 1:na], m[:r_sto][a] <= mg.storages[a].bounds.ub
        # Power bounds
        [h in 1:nh, s in 1:ns, a in 1:na], m[:p_dch][h,s,a]  <= mg.storages[a].α_p_dch * m[:r_sto][a]
        [h in 1:nh, s in 1:ns, a in 1:na], m[:p_ch][h,s,a] <= mg.storages[a].α_p_ch * m[:r_sto][a]
        # SoC bounds
        [h in 1:nh+1, s in 1:ns, a in 1:na], m[:soc][h,s,a]  <= mg.storages[a].α_soc_max * m[:r_sto][a]
        [h in 1:nh+1, s in 1:ns, a in 1:na], m[:soc][h,s,a]  >= mg.storages[a].α_soc_min * m[:r_sto][a]
        # State dynamics
        [h in 1:nh, s in 1:ns, a in 1:na], m[:soc][h+1,s,a]  == m[:soc][h,s,a] * (1. - mg.storages[a].η_self * mg.parameters.Δh) - (m[:p_dch][h,s,a] / mg.storages[a].η_dch - m[:p_ch][h,s,a] * mg.storages[a].η_ch) * mg.parameters.Δh
        # Initial and final states
        [s in 1:ns, a in 1:na], m[:soc][1,s,a]    == mg.storages[a].soc_ini * m[:r_sto][a]
        [s in 1:ns, a in 1:na], m[:soc][end,s,a]  >= m[:soc][1,s,a]
        end)
    end
    # Converters
    if !isempty(mg.converters)
        na = length(mg.converters)
        # Constraint
        @constraints(m, begin
        # Investment bounds
        [a in 1:na], m[:r_c][a] >= mg.converters[a].bounds.lb
        [a in 1:na], m[:r_c][a] <= mg.converters[a].bounds.ub
        # Power bounds
        [h in 1:nh, s in 1:ns, a in 1:na], m[:p_c][h,s,a]  <= m[:r_c][a]
        end)
    end
    # Grids
    if !isempty(mg.grids)
        na = length(mg.grids)
        # Constraint
        @constraints(m, begin
        # Power bounds
        [h in 1:nh, s in 1:ns, a in 1:na], m[:p_in][h,s,a]  <= mg.grids[a].powerMax
        [h in 1:nh, s in 1:ns, a in 1:na], m[:p_out][h,s,a] <= mg.grids[a].powerMax
        end)
    end
end
# Power balance
function add_power_balances!(m::Model, mg::Microgrid, ω::Scenarios)
    # !!! All the decision variables are defined positive !!!
    nh, ns = size(ω.demands[1].power,1), size(ω.demands[1].power,3)
    balance = (electricity = AffExpr.(zeros(nh,ns)), heat = AffExpr.(zeros(nh,ns)), hydrogen = AffExpr.(zeros(nh,ns)))
    # Demands
    for (k,a) in enumerate(mg.demands)
        if a.carrier isa Electricity
            add_to_expression!.(balance.electricity, ω.demands[k].power[:,1,:])
        elseif a.carrier isa Heat
            add_to_expression!.(balance.heat, ω.demands[k].power[:,1,:])
        elseif a.carrier isa Hydrogen
            add_to_expression!.(balance.hydrogen, ω.demands[k].power[:,1,:])
        end
    end
    # Generation
    for (k,a) in enumerate(mg.generations)
        if a.carrier isa Electricity
            add_to_expression!.(balance.electricity, .- m[:r_g][k] .* ω.generations[k].power[:,1,:])
        elseif a.carrier isa Heat
            add_to_expression!.(balance.heat, .- m[:r_g][k] .* ω.generations[k].power[:,1,:])
        elseif a.carrier isa Hydrogen
            add_to_expression!.(balance.hydrogen, .- m[:r_g][k] .* ω.generations[k].power[:,1,:])
        end
    end
    # Storages
    for (k,a) in enumerate(mg.storages)
        if a.carrier isa Electricity
            add_to_expression!.(balance.electricity, m[:p_ch][:,:,k] .- m[:p_dch][:,:,k])
        elseif a.carrier isa Heat
            add_to_expression!.(balance.heat,  m[:p_ch][:,:,k] .- m[:p_dch][:,:,k])
        elseif a.carrier isa Hydrogen
            add_to_expression!.(balance.hydrogen,  m[:p_ch][:,:,k] .- m[:p_dch][:,:,k])
        end
    end
    # Converters
    for (k,a) in enumerate(mg.converters)
        if a isa Heater
            add_to_expression!.(balance.electricity, m[:p_c][:,:,k])
            add_to_expression!.(balance.heat, .- m[:p_c][:,:,k] * a.η_E_H)
        elseif a isa Electrolyzer
            add_to_expression!.(balance.electricity, m[:p_c][:,:,k])
            add_to_expression!.(balance.heat, .- m[:p_c][:,:,k] * a.η_E_H)
            add_to_expression!.(balance.hydrogen, .- m[:p_c][:,:,k] * a.η_E_H2)
        elseif a isa FuelCell
            add_to_expression!.(balance.electricity, .- m[:p_c][:,:,k])
            add_to_expression!.(balance.heat, .- m[:p_c][:,:,k] * a.η_H2_H / a.η_H2_E)
            add_to_expression!.(balance.hydrogen, m[:p_c][:,:,k] / a.η_H2_E)
        end
    end
    # Grids
    for (k,a) in enumerate(mg.grids)
        if a.carrier isa Electricity
            add_to_expression!.(balance.electricity, .- m[:p_in][:,:,k] + m[:p_out][:,:,k])
        elseif a.carrier isa Heat
            add_to_expression!.(balance.heat, .- m[:p_in][:,:,k] + m[:p_out][:,:,k])
        elseif a.carrier isa Hydrogen
            add_to_expression!.(balance.hydrogen, .- m[:p_in][:,:,k] + m[:p_out][:,:,k])
        end
    end
    # Energy balance constraint
    @constraints(m, begin
        [h in 1:nh, s in 1:ns], balance.electricity[h,s] <= 0.
        [h in 1:nh, s in 1:ns], balance.heat[h,s] <= 0.
        [h in 1:nh, s in 1:ns], balance.hydrogen[h,s] == 0.
    end)
end
# Renewable share
function add_renewable_share!(m::Model, mg::Microgrid, ω::Scenarios, probabilities::Vector{Float64}, risk::AbstractRiskMeasure)
    nh, ns = size(ω.demands[1].power,1), size(ω.demands[1].power,3)
    total = zeros(ns)
    for (k,a) in enumerate(mg.demands)
        if a.carrier isa Electricity
            total .= total .+ sum(ω.demands[k].power[h,1,:] for h in 1:nh)
        elseif a.carrier isa Heat
            total .= total .+ sum(ω.demands[k].power[h,1,:] for h in 1:nh) ./ mg.converters[isin(mg.converters, Heater)[2]].η_E_H
        end
    end
    for (k,a) in enumerate(mg.grids)
        if a.carrier isa Electricity
            @expression(m, share[s in 1:ns], sum(m[:p_in][h,s,k] for h in 1:nh) - (1. - mg.parameters.renewable_share) * total[s])
        end
    end
    # Constraint according to CVaR
    @variables(m, begin
    ζ_s
    α_s[1:ns] >= 0.
    end)
    @constraints(m, begin
    [s in 1:ns], α_s[s] >= m[:share][s] - ζ_s
    ζ_s + 1 / (1 - beta(risk)) * sum(probabilities[s] * α_s[s] for s in 1:ns) <= 0.
    end)
end
# Objective
function add_objective!(m::Model, mg::Microgrid, ω::Scenarios, probabilities::Vector{Float64}, risk::AbstractRiskMeasure)
    ns = size(ω.demands[1].power,3)
    # CAPEX
    capex = compute_capex(m, mg, ω)
    # OPEX
    opex = compute_opex(m, mg, ω)
    # Objective according to the CVaR
    @variables(m, begin
    ζ_o
    α_o[1:ns] >= 0.
    end)
    @constraint(m, [s in 1:ns], α_o[s] >= capex + opex[s] - ζ_o)
    @objective(m, Min, ζ_o + 1 / (1 - beta(risk)) * sum(probabilities[s] * α_o[s] for s in 1:ns))
end
# Capex
function compute_capex(m::Model, mg::Microgrid, ω::Scenarios)
    cost = AffExpr(0.)
    # Generations
    for (k,a) in enumerate(mg.generations)
        add_to_expression!(cost, Γ(mg.parameters.τ, a.lifetime) * ω.generations[k].cost[1] * m[:r_g][k])
    end
    # Storages
    for (k,a) in enumerate(mg.storages)
        add_to_expression!(cost, Γ(mg.parameters.τ, a.lifetime) * ω.storages[k].cost[1] * m[:r_sto][k])
    end
    # Converters
    for (k,a) in enumerate(mg.converters)
        add_to_expression!(cost, Γ(mg.parameters.τ, a.lifetime) * ω.converters[k].cost[1] * m[:r_c][k])
    end
    return cost
end
# Grids
function compute_opex(m::Model, mg::Microgrid, ω::Scenarios)
    nh, ns = size(ω.demands[1].power,1), size(ω.demands[1].power,3)
    cost = AffExpr.(zeros(ns))
    for (k,a) in enumerate(mg.grids)
        add_to_expression!.(cost, sum((m[:p_in][h,:,k] .* ω.grids[k].cost_in[h,1,:] .- m[:p_out][h,:] .* ω.grids[k].cost_out[h,1,:]) .* mg.parameters.Δh  for h in 1:nh))
    end
    return cost
end


### Offline
function initialize_designer!(mg::Microgrid, designer::MILP, ω::Scenarios)
    # Preallocate
    preallocate!(mg, designer)

    # Scenario reduction from the optimization scenario pool
    if isa(designer.options.read_reduction, Nothing)
        println("Starting scenario reduction...")
        ω_reduced, probabilities = reduce(designer.options.reducer, ω)
        # Saving
        if !isa(designer.options.write_reduction, Nothing)
            save(designer.options.write_reduction, "scenarios", ω_reduced, "probabilities", probabilities)
        end
    else
        println("Reading scenario reduction from file...")
        ω_reduced = load(designer.options.read_reduction, "scenarios")
        probabilities = load(designer.options.read_reduction, "probabilities")
    end

    # Initialize model
    println("Building the model...")
    designer.model = build_model(mg, designer, ω_reduced, probabilities)

    # Compute investment decisions for the first year
    println("Starting optimization...")
    optimize!(designer.model)

    # Assign values
    for k in 1:length(mg.generations)
        designer.decisions.generations[k][1,:] .= value(designer.model[:r_g][k])
    end
    for k in 1:length(mg.storages)
        designer.decisions.storages[k][1,:] .= value(designer.model[:r_sto][k])
    end
    for k in 1:length(mg.converters)
        designer.decisions.converters[k][1,:] .= value(designer.model[:r_c][k])
    end

    # Save history
    designer.history = ω_reduced

     return designer
end

### Online
function compute_investment_decisions!(y::Int64, s::Int64, mg::Microgrid, designer::MILP)
    # TODO
end

### Utils
beta(risk::WorstCase) = 1. - 1e-6
beta(risk::Expectation) = 0.
beta(risk::CVaR) = risk.β
Γ(τ::Float64, lifetime::Union{Float64, Int64}) = τ * (τ + 1.) ^ lifetime / ((τ + 1.) ^ lifetime - 1.)
