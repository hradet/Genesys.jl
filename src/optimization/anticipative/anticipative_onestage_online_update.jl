#=
    Anticipative one-stage designer with anticipative controller and online
    updates
=#

mutable struct AnticipativeOneStageOnlineDesigner <: AbstractOneStageStochasticDesigner
    u::NamedTuple
    model::JuMP.Model
    parameters::Dict{String, Any}
    AnticipativeOneStageOnlineDesigner() = new()
end

#### Models ####
function onestage_milp_model(ld::Load, pv::Source, liion::Liion,
  designer::AnticipativeOneStageOnlineDesigner, grid::Grid, ω_optim::Scenarios,
  parameters::NamedTuple)

    # Parameters
    Γ_liion = (parameters.τ * (parameters.τ + 1.) ^ liion.lifetime) / ((parameters.τ + 1.) ^ liion.lifetime - 1.)
    Γ_pv = (parameters.τ * (parameters.τ + 1.) ^ pv.lifetime) / ((parameters.τ + 1.) ^ pv.lifetime - 1.)

    # Sets
    nh = size(ω_optim.values.ld_E,1) # Number of hours
    ns = size(ω_optim.values.ld_E,2) # Number of scenarios

    # Model definition
    m = Model(CPLEX.Optimizer)

    # Variables
    @variables(m, begin
    # Operation decision variables
     p_liion_ch[1:nh, 1:ns] <= 0.
     p_liion_dch[1:nh, 1:ns] >= 0.
     p_g_out[1:nh, 1:ns] <= 0.
     p_g_in[1:nh, 1:ns] >= 0.
     # Investment decision variables
     0 <= r_liion <= 1000
     0 <= r_pv <= 1000
     # Operation states variables
     soc_liion[1:nh+1, 1:ns]
    end)

    # Constraints
    @constraints(m, begin
    # Power bounds
    [h in 1:nh, s in 1:ns], p_liion_dch[h,s] <= liion.α_p_dch * r_liion
    [h in 1:nh, s in 1:ns], p_liion_ch[h,s] >= -liion.α_p_ch * r_liion
    [h in 1:nh, s in 1:ns], p_g_in[h,s] <= (1. - grid.τ_power) * maximum(ω_optim.values.ld_E)
    # SoC bounds
    [h in 1:nh+1, s in 1:ns], soc_liion[h,s] <= liion.α_soc_max * r_liion
    [h in 1:nh+1, s in 1:ns], soc_liion[h,s] >= liion.α_soc_min * r_liion
    # Investment bounds
    r_liion <= 1000.
    r_pv <= 1000.
    # State dynamic
    [h in 1:nh, s in 1:ns], soc_liion[h+1,s] == soc_liion[h,s] * (1 - liion.η_self * parameters.Δh) - (p_liion_ch[h,s] * liion.η_ch + p_liion_dch[h,s] / liion.η_dch) * parameters.Δh
    # Power balance
    [h in 1:nh, s in 1:ns], ω_optim.values.ld_E[h,s] - r_pv * ω_optim.values.pv_E[h,s] <= p_g_out[h,s] + p_g_in[h,s] + p_liion_ch[h,s] + p_liion_dch[h,s]
    # Self-sufficiency constraint
    self_constraint[s in 1:ns], sum(p_g_in[h,s] for h in 1:nh) <= (1. - grid.τ_energy) * sum(ω_optim.values.ld_E[h,s] for h in 1:nh)
    # Initial and final conditions
    [s in 1:ns], soc_liion[1,s] == liion.soc[1] * r_liion
    [s in 1:ns], soc_liion[nh,s] == soc_liion[1,s]
    end)

    # CAPEX
    capex = @expression(m, Γ_pv * ω_optim.values.C_pv[1] * r_pv + Γ_liion * ω_optim.values.C_liion[1] * r_liion)

    # OPEX
    if designer.parameters["risk"] == "esperance"
      #TODO : multiplier par proba du scenario au lieu de /ns...
      opex = @expression(m, sum((p_g_in[h,s] * ω_optim.values.C_grid_in[h,s] + p_g_out[h,s] * ω_optim.values.C_grid_out[h,s]) * parameters.Δh  for h in 1:nh, s in 1:ns) / ns)
    elseif designer.parameters["risk"]  == "cvar"
    else
      println("Unknown risk measure... Chose between 'esperance' or 'cvar' ")
    end

    # Objective
    @objective(m, Min, capex + opex)

    return m
end

#### Offline functions ####
function offline_optimization(ld::Load, pv::Source, liion::Liion,
      controller::AnticipativeController, designer::AnticipativeOneStageOnlineDesigner,
      grid::Grid, ω_optim::Scenarios, parameters::NamedTuple)
      # Parameters
      ny = size(ld.power_E,2) # number of simulation years
      ns = size(ld.power_E,3) # number of scenarios

      # Scenario reduction from the optimization scenario pool
      ω_anticipative = scenarios_reduction(designer, ω_optim)

      # Initialize model
      designer.model = controller.model = onestage_milp_model(ld, pv, liion, controller, designer, grid, ω_anticipative, parameters)

      # Compute both investment and operation decisions
      optimize!(designer.model)

      # Formatting variables to simulation
      # Operation decisions
      controller.u = (
      u_liion =  repeat(value.(controller.model[:p_liion_ch] + controller.model[:p_liion_dch]), 1, 1, ns),
      )
      # Investment controls
      designer.u = (
      u_liion = repeat(vcat(value.(designer.model[:r_liion]), zeros(ny-1,1)),1, ns),
      u_pv = repeat(vcat(value.(designer.model[:r_pv]), zeros(ny-1,1)),1, ns),
      )
end

#### Online functions ####
function compute_investment_decisions(y::Int64, s::Int64, ld::Load, pv::Source,
    liion::Liion, grid::Grid, designer::AnticipativeOneStageOnlineDesigner, ω_optim::Scenarios, parameters::NamedTuple)
    # Parameters
    ϵ = 0.1
    ny = size(ld.power_E,2)
    ns = size(ld.power_E,3)

    if liion.soh[end,y,s] < ϵ
        # Update reduction range
        designer.parameters["idx_years"] = y+1:ny
        # Scenario reduction from the scenario pool
        ω_anticipative = scenarios_reduction(designer, ω_optim)
        # Initialize model with the new values
        designer.model = controller.model = onestage_milp_model(ld, pv, liion, controller, designer, grid, ω_anticipative, parameters)
        # Compute both investment and operation decisions
        # Fix PV design
        fix.(designer.model[:r_pv], pv.powerMax[y], force = true)
        # Optimize
        optimize!(designer.model)
        # Formatting variables to simulation
        # Operation decisions
        controller.u.u_liion[:,y+1:end,s] = repeat(value.(controller.model[:p_liion_ch] + controller.model[:p_liion_dch]), 1, ceil(Int64,(ny-y)/size(controller.model[:p_liion_ch],2)))
        # Investment controls
        designer.u.u_liion[y,s] = value.(designer.model[:r_liion])
    end
end
