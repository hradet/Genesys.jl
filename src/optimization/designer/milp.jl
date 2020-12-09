#=
    Designer based on the equivalent annual cost (EAC) with multiple scenarios
=#

mutable struct MILPOptions
  solver::Module
  risk_measure::String
  q::Float64
  reducer::AbstractScenariosReducer
  share_constraint::String
  reopt::Bool

  MILPOptions(; solver = CPLEX,
                risk_measure = "expectation",
                q = 0.95,
                reducer = KmeansReducer(),
                share_constraint = "hard",
                reopt=false) =
                new(solver, risk_measure, q, reducer, share_constraint, reopt)
end

mutable struct MILP <: AbstractDesigner
    options::MILPOptions
    u::NamedTuple
    model::JuMP.Model
    history::AbstractScenarios

    MILP(; options = MILPOptions()) = new(options)
end

### Models
function build_model(des::DistributedEnergySystem, designer::MILP, ω::Scenarios, probabilities::Vector{Float64})
    # Sets
    nh = size(ω.ld_E.power,1) # Number of hours
    ns = size(ω.ld_E.power,3) # Number of scenarios

    # Initialize
    isa(des.ld_E, Load) ? ld_E = ω.ld_E.power : ld_E = zeros(nh, 1, ns)
    isa(des.ld_H, Load) ? ld_H = ω.ld_H.power : ld_H = zeros(nh, 1, ns)
    isa(des.liion, Liion) ? liion = des.liion : liion = Liion()
    isa(des.tes, ThermalSto) ? tes = des.tes : tes = ThermalSto()
    isa(des.h2tank, H2Tank) ? h2tank = des.h2tank : h2tank = H2Tank()
    isa(des.elyz, Electrolyzer) ? elyz = des.elyz : elyz = Electrolyzer()
    isa(des.fc, FuelCell) ? fc = des.fc : fc = FuelCell()
    isa(des.heater, Heater) ? heater = des.heater : heater = Heater()
    isa(des.pv, Source) ? pv = des.pv : pv = Source()
    isa(des.grid, Grid) ? grid = des.grid : grid = Grid()

    # Model definition
    m = Model(designer.options.solver.Optimizer)

    # Build model
    # Variables
    @variables(m, begin
    # Operation decisions variables
    p_liion_ch[1:nh, 1:ns] <= 0.
    p_liion_dch[1:nh, 1:ns] >= 0.
    p_tes_ch[1:nh, 1:ns] <= 0.
    p_tes_dch[1:nh, 1:ns] >= 0.
    p_h2tank_ch[1:nh, 1:ns] <= 0.
    p_h2tank_dch[1:nh, 1:ns] >= 0.
    p_elyz_E[1:nh, 1:ns] <= 0.
    p_fc_E[1:nh, 1:ns] >= 0.
    p_heater_E[1:nh, 1:ns] <= 0.
    p_g_out[1:nh, 1:ns] <= 0.
    p_g_in[1:nh, 1:ns] >= 0.
    # Investment decisions variables
    0 <= r_liion <= (isa(des.liion, Liion) ? 1000 : 0)
    0 <= r_tes <= (isa(des.tes, ThermalSto) ? 1000 : 0)
    0 <= r_h2tank <= (isa(des.h2tank, H2Tank) ? 50000 : 0)
    0 <= r_elyz <= (isa(des.elyz, Electrolyzer) ? 50 : 0)
    0 <= r_fc <= (isa(des.fc, FuelCell) ? 50 : 0)
    0 <= r_pv <= (isa(des.pv, Source) ? 1000 : 0)
    # Operation state variables
    soc_liion[1:nh+1, 1:ns]
    soc_tes[1:nh+1, 1:ns]
    soc_h2tank[1:nh+1, 1:ns]
    end)
    # Constraints
    @constraints(m, begin
    # Power bounds
    [h in 1:nh, s in 1:ns], p_liion_dch[h,s] <= liion.α_p_dch * r_liion
    [h in 1:nh, s in 1:ns], p_liion_ch[h,s] >= -liion.α_p_ch * r_liion
    [h in 1:nh, s in 1:ns], p_tes_dch[h,s] <= tes.α_p_dch * r_tes
    [h in 1:nh, s in 1:ns], p_tes_ch[h,s] >= -tes.α_p_ch * r_tes
    [h in 1:nh, s in 1:ns], p_h2tank_dch[h,s] <= h2tank.α_p_dch * r_h2tank
    [h in 1:nh, s in 1:ns], p_h2tank_ch[h,s] >= -h2tank.α_p_ch * r_h2tank
    [h in 1:nh, s in 1:ns], p_elyz_E[h,s] >= -r_elyz
    [h in 1:nh, s in 1:ns], p_fc_E[h,s] <= r_fc
    [h in 1:nh, s in 1:ns], p_heater_E[h,s] >= -heater.powerMax_ini
    [h in 1:nh, s in 1:ns], p_g_in[h,s] <= grid.powerMax
    [h in 1:nh, s in 1:ns], p_g_out[h,s] >= -grid.powerMax
    # SoC bounds
    [h in 1:nh+1, s in 1:ns], soc_liion[h,s] <= liion.α_soc_max * r_liion
    [h in 1:nh+1, s in 1:ns], soc_liion[h,s] >= liion.α_soc_min * r_liion
    [h in 1:nh+1, s in 1:ns], soc_tes[h,s] <= tes.α_soc_max * r_tes
    [h in 1:nh+1, s in 1:ns], soc_tes[h,s] >= tes.α_soc_min * r_tes
    [h in 1:nh+1, s in 1:ns], soc_h2tank[h,s] <= h2tank.α_soc_max * r_h2tank
    [h in 1:nh+1, s in 1:ns], soc_h2tank[h,s] >= h2tank.α_soc_min * r_h2tank
    # State dynamics
    [h in 1:nh, s in 1:ns], soc_liion[h+1,s] == soc_liion[h,s] * (1. - liion.η_self * des.parameters.Δh) - (p_liion_ch[h,s] * liion.η_ch + p_liion_dch[h,s] / liion.η_dch) * des.parameters.Δh
    [h in 1:nh, s in 1:ns], soc_tes[h+1,s] == soc_tes[h,s] * (1. - tes.η_self * des.parameters.Δh) - (p_tes_ch[h,s] * tes.η_ch + p_tes_dch[h,s] / tes.η_dch) * des.parameters.Δh
    [h in 1:nh, s in 1:ns], soc_h2tank[h+1,s] == soc_h2tank[h,s] * (1. - h2tank.η_self * des.parameters.Δh) - (p_h2tank_ch[h,s] * h2tank.η_ch + p_h2tank_dch[h,s] / h2tank.η_dch) * des.parameters.Δh
    # Initial and final states
    [s in 1:ns], soc_liion[1,s] == liion.soc_ini * r_liion
    [s in 1:ns], soc_liion[end,s] >= soc_liion[1,s]
    [s in 1:ns], soc_tes[1,s] == tes.soc_ini * r_tes
    [s in 1:ns], soc_tes[end,s] >= soc_tes[1,s]
    [s in 1:ns], soc_h2tank[1,s] == h2tank.soc_ini * r_h2tank
    [s in 1:ns], soc_h2tank[end,s] >= soc_h2tank[1,s]
    # Power balances
    [h in 1:nh, s in 1:ns], ld_E[h,1,s] <= r_pv * ω.pv.power[h,1,s] + p_liion_ch[h,s] + p_liion_dch[h,s] + p_elyz_E[h,s] + p_fc_E[h,s] + p_heater_E[h,s] + p_g_in[h,s] + p_g_out[h,s]
    [h in 1:nh, s in 1:ns], ld_H[h,1,s] <= p_tes_ch[h,s]  + p_tes_dch[h,s] - elyz.η_E_H * p_elyz_E[h,s] + fc.η_H2_H / fc.η_H2_E * p_fc_E[h,s] - heater.η_E_H * p_heater_E[h,s]
    [h in 1:nh, s in 1:ns], 0. == p_h2tank_ch[h,s] + p_h2tank_dch[h,s] - elyz.η_E_H2 * p_elyz_E[h,s] - p_fc_E[h,s] / fc.η_H2_E
    end)

    # Share of renewables constraint
    if designer.options.share_constraint == "hard"
        @constraint(m, [s in 1:ns], sum(p_g_in[h,s] - (1. - des.parameters.τ_share) * (ld_E[h,1,s] + ld_H[h,1,s] / heater.η_E_H) for h in 1:nh) <= 0.)
    elseif designer.options.share_constraint == "soft"
        @constraint(m, sum(probabilities[s] * (p_g_in[h,s] - (1. - des.parameters.τ_share) * (ld_E[h,1,s] + ld_H[h,1,s] / heater.η_E_H)) for h in 1:nh, s in 1:ns) <= 0.)
    end

    # CAPEX
    # Annualized factor
    Γ_liion = (des.parameters.τ * (des.parameters.τ + 1.) ^ liion.lifetime) / ((des.parameters.τ + 1.) ^ liion.lifetime - 1.)
    Γ_tes = (des.parameters.τ * (des.parameters.τ + 1.) ^ tes.lifetime) / ((des.parameters.τ + 1.) ^ tes.lifetime - 1.)
    Γ_h2tank = (des.parameters.τ * (des.parameters.τ + 1.) ^ h2tank.lifetime) / ((des.parameters.τ + 1.) ^ h2tank.lifetime - 1.)
    Γ_elyz = (des.parameters.τ * (des.parameters.τ + 1.) ^ elyz.lifetime) / ((des.parameters.τ + 1.) ^ elyz.lifetime - 1.)
    Γ_fc = (des.parameters.τ * (des.parameters.τ + 1.) ^ fc.lifetime) / ((des.parameters.τ + 1.) ^ fc.lifetime - 1.)
    Γ_pv = (des.parameters.τ * (des.parameters.τ + 1.) ^ pv.lifetime) / ((des.parameters.τ + 1.) ^ pv.lifetime - 1.)
    @expression(m, capex, Γ_liion * ω.liion.cost[1] * r_liion + Γ_tes * ω.tes.cost[1] * r_tes + Γ_h2tank * ω.h2tank.cost[1] * r_h2tank + Γ_elyz * ω.elyz.cost[1] * r_elyz + Γ_fc * ω.fc.cost[1] * r_fc + Γ_pv * ω.pv.cost[1] * r_pv)

    # OPEX
    @expression(m, opex[s in 1:ns], sum((p_g_in[h,s] * ω.grid.cost_in[h,1,s] + p_g_out[h,s] * ω.grid.cost_out[h,1,s]) * des.parameters.Δh  for h in 1:nh))

    # Objective
    if designer.options.risk_measure == "expectation"
        @objective(m, Min, capex + sum(probabilities[s] * opex[s] for s in 1:ns))
    elseif designer.options.risk_measure == "cvar"
        @variable(m, ζ)
        @variable(m, α[1:ns] >= 0.)
        @constraint(m, [s in 1:ns], α[s] >= capex + opex[s] - ζ)
        @objective(m, Min, ζ + 1 / (1 - designer.options.q) * sum(probabilities[s] * α[s] for s in 1:ns))
    else
      println("Unknown risk measure... Chose between 'esperance' or 'cvar' ")
    end

    return m
end

### Offline
function initialize_designer!(des::DistributedEnergySystem, designer::MILP, ω::Scenarios{Array{DateTime,3}, Array{Float64,3}, Array{Float64,2}})
    # Preallocate
    preallocate!(designer, des.parameters.ny, des.parameters.ns)

    # Scenario reduction from the optimization scenario pool
    println("Starting scenario reduction...")
    ω_reduced, probabilities = reduce(designer.options.reducer, ω)

    # Initialize model
    println("Building the model...")
    designer.model = build_model(des, designer, ω_reduced, probabilities)

    # Compute investment decisions for the first year
    println("Starting optimization...")
    optimize!(designer.model)

    # Assign values
    designer.u.pv[1,:] .= value(designer.model[:r_pv])
    designer.u.liion[1,:] .= value(designer.model[:r_liion])
    designer.u.h2tank[1,:] .= value(designer.model[:r_h2tank])
    designer.u.elyz[1,:] .= value(designer.model[:r_elyz])
    designer.u.fc[1,:] .= value(designer.model[:r_fc])
    designer.u.tes[1,:] .= value(designer.model[:r_tes])

    # Save history
    designer.history = ω_reduced

     return designer
end

### Online
function compute_investment_decisions!(y::Int64, s::Int64, des::DistributedEnergySystem, designer::MILP)
    ϵ = 0.1

    # TODO : fix reoptimize function !!
    if designer.options.reopt && y != 1
        # Do we need to reoptimize ?
        (isa(des.liion, Liion) && des.liion.soh[end,y,s] < ϵ) || (isa(des.elyz, Electrolyzer) && des.elyz.soh[end,y,s] < ϵ) || (isa(des.fc, FuelCell) && des.fc.soh[end,y,s] < ϵ) ? nothing : return

        # Scenario reduction from the optimization scenario pool
        ω = scenarios_reduction(designer.history, 1:des.parameters.nh, y:des.parameters.ny, 1)

        # Build model
        designer.model = build_model(des, designer, ω)

        # Fix variables
        isa(des.pv, Source) ? fix(designer.model[:r_pv], des.pv.powerMax[y,s], force=true) : nothing
        isa(des.h2tank, H2Tank) ? fix(designer.model[:r_h2tank], des.h2tank.Erated[y,s], force=true) : nothing
        isa(des.tes, ThermalSto) ? fix(designer.model[:r_tes], des.tes.Erated[y,s], force=true) : nothing
        isa(des.liion, Liion) && des.liion.soh[end,y,s] >= ϵ ? fix(designer.model[:r_liion], des.liion.Erated[y,s], force=true) : nothing
        isa(des.elyz, Electrolyzer) && des.elyz.soh[end,y,s] >= ϵ ? fix(designer.model[:r_elyz], des.elyz.powerMax[y,s], force=true) : nothing
        isa(des.fc, FuelCell) && des.fc.soh[end,y,s] >= ϵ ? fix(designer.model[:r_fc], des.fc.powerMax[y,s], force=true) : nothing

        # Compute investment decisions
        optimize!(designer.model)

        # Preallocate and assigned values
        isa(des.liion, Liion) && des.liion.soh[end,y,s] < ϵ ? designer.u.liion[y,s] = value(designer.model[:r_liion]) : nothing
        isa(des.elyz, Electrolyzer) && des.elyz.soh[end,y,s] < ϵ ? designer.u.elyz[y,s] = value(designer.model[:r_elyz]) : nothing
        isa(des.fc, FuelCell) && des.fc.soh[end,y,s] < ϵ ? designer.u.fc[y,s] = value(designer.model[:r_fc]) : nothing

    else
        isa(des.liion, Liion) && des.liion.soh[end,y,s] < ϵ ? designer.u.liion[y,s] = designer.u.liion[1,s] : nothing
        isa(des.elyz, Electrolyzer) && des.elyz.soh[end,y,s] < ϵ ? designer.u.elyz[y,s] = designer.u.elyz[1,s] : nothing
        isa(des.fc, FuelCell) && des.fc.soh[end,y,s] < ϵ ? designer.u.fc[y,s] = designer.u.fc[1,s] : nothing
    end
end
