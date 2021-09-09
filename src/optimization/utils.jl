# Risk measures
abstract type AbstractRiskMeasure end
# Expectation
struct Expectation <: AbstractRiskMeasure end
# Worst case
struct WorstCase <: AbstractRiskMeasure end
# Conditional value at risk
struct CVaR <: AbstractRiskMeasure
    β::Float64

    function CVaR(β::Float64)
        @assert 0. <= β <= 1. "β must be in [0,1]"
        return new(β)
    end
end

conditional_value_at_risk(support::Array{Float64,1}, probabilities::Array{Float64,1}, risk::WorstCase) = conditional_value_at_risk(support, probabilities, 0.)
conditional_value_at_risk(support::Array{Float64,1}, probabilities::Array{Float64,1}, risk::Expectation) = conditional_value_at_risk(support, probabilities, 1.)
conditional_value_at_risk(support::Array{Float64,1}, probabilities::Array{Float64,1}, risk::CVaR) = conditional_value_at_risk(support, probabilities, 1. - risk.β)

# From https://github.com/jaantollander/ConditionalValueAtRisk
function conditional_value_at_risk(support::Array{Float64,1}, probabilities::Array{Float64,1}, α::Float64)
    # Value at risk
    var = value_at_risk(support, probabilities, α)
    # Tail values
    if α == 0.
        return var
    else
        tail = support .< var
        return (sum(probabilities[tail] .* support[tail]) - (sum(probabilities[tail]) - α) * var) / α
    end
end

function value_at_risk(support::Array{Float64,1}, probabilities::Array{Float64,1}, α::Float64)
    i = findfirst(cumsum(probabilities[sortperm(support)]) .>= α)
    if i == nothing
        return sort(support)[end]
    else
        return sort(support)[i]
    end
end

# MILP functions
# Decisions
function add_investment_decisions!(m::Model, generations::Vector{AbstractGeneration})
    if !isempty(generations)
        na = length(generations)
        @variable(m, r_g[1:na])
    end
end
function add_investment_decisions!(m::Model, storages::Vector{AbstractStorage})
    if !isempty(storages)
        na = length(storages)
        @variable(m, r_sto[1:na])
    end
end
function add_investment_decisions!(m::Model, converters::Vector{AbstractConverter})
    if !isempty(converters)
        na = length(converters)
        @variable(m, r_c[1:na])
    end
end
function fix_investment_decisions!(m::Model, generations::Vector{Float64}, storages::Vector{Float64}, converters::Vector{Float64})
    # Generation
    if !isempty(generations)
        fix.(m[:r_g], generations)
    end
    # Storages
    if !isempty(storages)
        fix.(m[:r_sto], storages)
    end
    # Converters
    if !isempty(converters)
        fix.(m[:r_c], converters)
    end
end
function add_operation_decisions!(m::Model, demands::Vector{AbstractDemand}, nh::Int64, ns::Int64)
    if !isempty(demands)
        na = length(demands)
        @variables(m, begin
        p_d[1:nh, 1:ns, 1:na]
        end)
    end
end
function add_operation_decisions!(m::Model, generations::Vector{AbstractGeneration}, nh::Int64, ns::Int64)
    if !isempty(generations)
        na = length(generations)
        @variables(m, begin
        p_g[1:nh, 1:ns, 1:na]
        end)
    end
end
function add_operation_decisions!(m::Model, storages::Vector{AbstractStorage}, nh::Int64, ns::Int64)
    if !isempty(storages)
        na = length(storages)
        @variables(m, begin
        p_ch[1:nh, 1:ns, 1:na]   >= 0.
        p_dch[1:nh, 1:ns, 1:na]  >= 0.
        soc[1:nh+1, 1:ns, 1:na]
        end)
    end
end
function add_operation_decisions!(m::Model, converters::Vector{AbstractConverter}, nh::Int64, ns::Int64)
    if !isempty(converters)
        na = length(converters)
        @variable(m, p_c[1:nh, 1:ns, 1:na] >= 0.)
    end
end
function add_operation_decisions!(m::Model, grids::Vector{AbstractGrid}, nh::Int64, ns::Int64)
    if !isempty(grids)
        na = length(grids)
        @variables(m, begin
        p_in[1:nh, 1:ns, 1:na]   >= 0.
        p_out[1:nh, 1:ns, 1:na]  >= 0.
        end)
    end
end
# Investment bounds
function add_investment_constraints!(m::Model, generations::Vector{AbstractGeneration})
    if !isempty(generations)
        na = length(generations)
        @constraints(m, begin
        [a in 1:na], m[:r_g][a] >= generations[a].bounds.lb
        [a in 1:na], m[:r_g][a] <= generations[a].bounds.ub
        end)
    end
end
function add_investment_constraints!(m::Model, storages::Vector{AbstractStorage})
    if !isempty(storages)
        na = length(storages)
        @constraints(m, begin
        [a in 1:na], m[:r_sto][a] >= storages[a].bounds.lb
        [a in 1:na], m[:r_sto][a] <= storages[a].bounds.ub
        end)
    end
end
function add_investment_constraints!(m::Model, converters::Vector{AbstractConverter})
    if !isempty(converters)
        na = length(converters)
        @constraints(m, begin
        [a in 1:na], m[:r_c][a] >= converters[a].bounds.lb
        [a in 1:na], m[:r_c][a] <= converters[a].bounds.ub
        end)
    end
end
# Technical constraint
function add_technical_constraints!(m::Model, storages::Vector{AbstractStorage}, Δh::Int64, nh::Int64, ns::Int64)
    if !isempty(storages)
        na = length(storages)
        @constraints(m, begin
        # Power bounds
        [h in 1:nh, s in 1:ns, a in 1:na], m[:p_dch][h,s,a] <= storages[a].α_p_dch * m[:r_sto][a]
        [h in 1:nh, s in 1:ns, a in 1:na], m[:p_ch][h,s,a]  <= storages[a].α_p_ch * m[:r_sto][a]
        # SoC bounds
        [h in 1:nh+1, s in 1:ns, a in 1:na], m[:soc][h,s,a] <= storages[a].α_soc_max * m[:r_sto][a]
        [h in 1:nh+1, s in 1:ns, a in 1:na], m[:soc][h,s,a] >= storages[a].α_soc_min * m[:r_sto][a]
        # State dynamics
        [h in 1:nh, s in 1:ns, a in 1:na], m[:soc][h+1,s,a] == m[:soc][h,s,a] * (1. - storages[a].η_self * Δh) - (m[:p_dch][h,s,a] / storages[a].η_dch - m[:p_ch][h,s,a] * storages[a].η_ch) * Δh
        # Initial and final states
        [s in 1:ns, a in 1:na], m[:soc][1,s,a] == storages[a].soc_ini * m[:r_sto][a]
        end)
    end
end
function add_technical_constraints!(m::Model, converters::Vector{AbstractConverter}, nh::Int64, ns::Int64)
    if !isempty(converters)
        na = length(converters)
        @constraints(m, begin
        # Power bounds
        [h in 1:nh, s in 1:ns, a in 1:na], m[:p_c][h,s,a]  <= m[:r_c][a]
        end)
    end
end
function add_technical_constraints!(m::Model, grids::Vector{AbstractGrid}, nh::Int64, ns::Int64)
    if !isempty(grids)
        na = length(grids)
        @constraints(m, begin
        # Power bounds
        [h in 1:nh, s in 1:ns, a in 1:na], m[:p_in][h,s,a]  <= grids[a].powerMax
        [h in 1:nh, s in 1:ns, a in 1:na], m[:p_out][h,s,a] <= grids[a].powerMax
        end)
    end
end
# Periodicity constraint
function add_periodicity_constraints!(m::Model, storages::Vector{AbstractStorage}, ns::Int64)
    # Storages
    if !isempty(storages)
        na = length(storages)
        @constraints(m, begin
        # Final states
        [s in 1:ns, a in 1:na], m[:soc][end,s,a]  >= m[:soc][1,s,a]
        end)
    end
end
# Power balance
function add_power_balance!(m::Model, mg::Microgrid, ω::Scenarios, type::DataType, nh::Int64, ns::Int64; ispnet::Bool=false)
    # !!! All the decision variables are defined positive !!!
    balance = AffExpr.(zeros(nh,ns))
    # Demands and generation
    if !ispnet
        for (k,a) in enumerate(mg.demands)
            if a.carrier isa type
                add_to_expression!.(balance, ω.demands[k].power[:,1,:])
            end
        end
        # Generation
        for (k,a) in enumerate(mg.generations)
            if a.carrier isa type
                add_to_expression!.(balance, .- m[:r_g][k] .* ω.generations[k].power[:,1,:])
            end
        end
    else
        for (k,a) in enumerate(mg.demands)
            if a.carrier isa type
                add_to_expression!.(balance, m[:p_d][:,:,k])
            end
        end
        # Generation
        for (k,a) in enumerate(mg.generations)
            if a.carrier isa type
                add_to_expression!.(balance, .- m[:p_g][:,:,k])
            end
        end
    end
    # Storages
    for (k,a) in enumerate(mg.storages)
        if a.carrier isa type
            add_to_expression!.(balance, m[:p_ch][:,:,k] .- m[:p_dch][:,:,k])
        end
    end
    # Converters
    for (k,a) in enumerate(mg.converters)
        if type == Electricity
            if a isa Heater
                add_to_expression!.(balance, m[:p_c][:,:,k])
            elseif a isa Electrolyzer
                add_to_expression!.(balance, m[:p_c][:,:,k])
            elseif a isa FuelCell
                add_to_expression!.(balance, .- m[:p_c][:,:,k])
            end
        elseif type == Heat
            if a isa Heater
                add_to_expression!.(balance, .- m[:p_c][:,:,k] * a.η_E_H)
            elseif a isa Electrolyzer
                add_to_expression!.(balance, .- m[:p_c][:,:,k] * a.η_E_H)
            elseif a isa FuelCell
                add_to_expression!.(balance, .- m[:p_c][:,:,k] * a.η_H2_H / a.η_H2_E)
            end
        elseif type == Hydrogen
            if a isa Electrolyzer
                add_to_expression!.(balance, .- m[:p_c][:,:,k] * a.η_E_H2)
            elseif a isa FuelCell
                add_to_expression!.(balance, m[:p_c][:,:,k] / a.η_H2_E)
            end
        end
    end
    # Grids
    for (k,a) in enumerate(mg.grids)
        if a.carrier isa type
            add_to_expression!.(balance, .- m[:p_in][:,:,k] + m[:p_out][:,:,k])
        end
    end
    # Energy balance constraint
    if type == Electricity
        @constraint(m, balance .<= 0.)
    elseif type == Heat
        @constraint(m, balance .<= 0.)
    elseif type == Hydrogen
        @constraint(m, balance .== 0.)
    end
end
# Renewable share
function add_renewable_share!(m::Model, mg::Microgrid, ω::Scenarios, probabilities::Vector{Float64}, risk::AbstractRiskMeasure, nh::Int64, ns::Int64)
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
function add_design_objective!(m::Model, mg::Microgrid, ω::Scenarios, probabilities::Vector{Float64}, risk::AbstractRiskMeasure, nh::Int64, ns::Int64)
    # CAPEX
    capex = compute_capex(m, mg, ω)
    # OPEX
    opex = compute_opex(m, mg, ω, nh, ns)
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
function compute_opex(m::Model, mg::Microgrid, ω::Scenarios, nh::Int64, ns::Int64)
    cost = AffExpr.(zeros(ns))
    for (k,a) in enumerate(mg.grids)
        add_to_expression!.(cost, sum((m[:p_in][h,:,k] .* ω.grids[k].cost_in[h,1,:] .- m[:p_out][h,:] .* ω.grids[k].cost_out[h,1,:]) .* mg.parameters.Δh  for h in 1:nh))
    end
    return cost
end
