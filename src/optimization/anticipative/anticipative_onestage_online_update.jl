#=
    Anticipative one-stage designer with anticipative controller and online
    updates
=#

mutable struct AnticipativeOneStageOnlineDesigner <: AbstractDesigner
    u::NamedTuple
    horizon::Int64
    model::JuMP.Model
    AnticipativeOneStageOnlineDesigner() = new()
end

#### Models ####
function onestage_milp_model(ld::Load, pv::Source, liion::Liion,
  controller::AnticipativeController, designer::AnticipativeOneStageOnlineDesigner,
  grid::Grid, ω_optim::Scenarios, parameters::NamedTuple)

    # Parameters
    crf_liion = (parameters.τ * (parameters.τ + 1) ^ liion.lifetime) / ((parameters.τ + 1) ^ liion.lifetime - 1)
    crf_pv = (parameters.τ * (parameters.τ + 1) ^ pv.lifetime) / ((parameters.τ + 1) ^ pv.lifetime - 1)

    # Sets
    nh = length(1:parameters.Δh:controller.horizon) # Number of hours
    ns = size(ω_optim.ld_E,2) # Number of one-year scenarios

    # Model definition
    m = Model(CPLEX.Optimizer)

    # Variables
    # Operation control variables
    @variables(m, begin
    # Liion
     p_liion_ch[1:nh,1:ns] <= 0.
     p_liion_dch[1:nh,1:ns] >= 0.
     # Recourse
     p_g_out[1:nh,1:ns] <= 0.
     p_g_in[1:nh,1:ns] >= 0.
    end)

    # Investment control variables
    @variables(m, begin
    r_liion >= 0.
    r_pv >= 0.
    end)

    # Operation state variables
    @variables(m, begin
    soc_liion[1:nh+1,1:ns]
    end)

    # Operation constraints bounds
    @constraints(m, begin
    # Controls
    # Liion
    [h in 1:nh, s in 1:ns], p_liion_dch[h,s] <= liion.α_p_dch * r_liion
    [h in 1:nh, s in 1:ns], p_liion_ch[h,s] >= -liion.α_p_ch * r_liion

    # State
    # Liion
    [h in 1:nh+1, s in 1:ns], soc_liion[h,s] <= liion.α_soc_max * r_liion
    [h in 1:nh+1, s in 1:ns], soc_liion[h,s] >= liion.α_soc_min * r_liion
    end)

    # Investment constraints bounds
    @constraints(m, begin
    # Controls
    r_liion <= 1000.
    r_pv <= 1000.
    end)

    # Operation constraints dynamics and recourse
    @constraints(m, begin
    # Dynamics
    # Liion
    [h in 1:nh, s in 1:ns], soc_liion[h+1,s] == soc_liion[h,s] * (1 - liion.η_self * parameters.Δh) - (p_liion_ch[h,s] * liion.η_ch + p_liion_dch[h,s] / liion.η_dch) * parameters.Δh

    # Recourse
    [h in 1:nh, s in 1:ns], ω_optim.ld_E[h,s] - r_pv * ω_optim.pv_E[h,s] <= p_g_out[h,s] + p_g_in[h,s] + p_liion_ch[h,s] + p_liion_dch[h,s]
    end)

    # Initial and final conditions
    @constraints(m, begin
    [s in 1:ns], soc_liion[1,s] == liion.soc[1] * r_liion
    [s in 1:ns], soc_liion[end,s] == soc_liion[1,s]
    end)

    # Grid constraints
    @constraints(m, begin
    power_constraint, p_g_in .<= (1. - grid.τ_power) * maximum(ω_optim.ld_E)
    self_constraint[s in 1:ns], sum(p_g_in[h,s] for h in 1:nh) <= (1. - grid.τ_energy) * sum(ω_optim.ld_E[h,s] for h in 1:nh)
    end)

    # CAPEX
    capex = @expression(m, crf_pv * ω_optim.C_pv[1] * r_pv + crf_liion * ω_optim.C_liion[1] * r_liion)

    # OPEX
    opex = @expression(m, sum((p_g_in[h,s] * ω_optim.C_grid_in[h,s] + p_g_out[h,s] * ω_optim.C_grid_out[h,s]) * parameters.Δh  for h in 1:nh, s in 1:ns) / ns)

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

      # Selection of one-year scenarios from the optimization dataset
      ω_milp = Scenarios(ω_optim.timestamp,ω_optim.ld_E[:,:,1],nothing,ω_optim.pv_E[:,:,1],
      ω_optim.C_pv[:,1],ω_optim.C_liion[:,1],nothing,nothing,nothing,nothing,nothing,
      ω_optim.C_grid_in[:,:,1],ω_optim.C_grid_out[:,:,1])

      # Initialize model
      designer.model = controller.model = onestage_milp_model(
      ld, pv, liion, controller, designer, grid, ω_milp, parameters)

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
    liion::Liion, grid::Grid, controller::AnticipativeController,
    designer::AnticipativeOneStageOnlineDesigner, ω_optim::Scenarios, parameters::NamedTuple)
    # Parameters
    ϵ = 0.1
    ny = size(ld.power_E,2)
    ns = size(ld.power_E,3)

    if liion.soh[end,y,s] < ϵ
        # Compute forecast scenarios
        ω_forecast = Scenarios(ω_optim.timestamp,ω_optim.ld_E[:,y+1:end,s],nothing,ω_optim.pv_E[:,y+1:end,s],
        ω_optim.C_pv[y,s],ω_optim.C_liion[y,s],nothing,nothing,nothing,nothing,nothing,
        ω_optim.C_grid_in[:,y+1:end,s],ω_optim.C_grid_out[:,y+1:end,s])
        # Initialize model
        designer.model = controller.model = onestage_milp_model(
        ld, pv, liion, controller, designer, grid, ω_forecast, parameters)
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
