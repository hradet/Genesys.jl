#=
    Model predictive control controller
=#

mutable struct MPCController <: AbstractController
    u::NamedTuple
    horizon::Int64
    markovchains::NamedTuple
    model::JuMP.Model
    MPCController() = new()
end

#### Models ####
# Simple
function mpc_model(ld::Load, pv::Source, liion::Liion, controller::MPCController,
     grid::Grid, parameters::NamedTuple)

     # Sets
     nh = length(1:parameters.Δh:controller.horizon) # Number of hours

     # Model definition
     m = Model(CPLEX.Optimizer)
     set_optimizer_attribute(m,"CPX_PARAM_SCRIND", 0)

     # Variables
     @variables(m, begin
     # Operation decisions variables
     p_liion_ch[1:nh] <= 0.
     p_liion_dch[1:nh] >= 0.
     p_g_out[1:nh] <= 0.
     p_g_in[1:nh] >= 0.
     # State variables
     soc_liion[1:nh+1]
     # Variables to be updated
     p_net[1:nh]
     E_liion
     soc_init
     end)

     # Constraints
     @constraints(m, begin
     # Power bounds
     [h in 1:nh], p_liion_dch[h] <= liion.α_p_dch * E_liion
     [h in 1:nh], p_liion_ch[h] >= -liion.α_p_ch * E_liion
     # SoC bounds
     [h in 1:nh+1], soc_liion[h] <= liion.α_soc_max * E_liion
     [h in 1:nh+1], soc_liion[h] >= liion.α_soc_min * E_liion
     # Dynamics
     [h in 1:nh], soc_liion[h+1] == soc_liion[h] * (1 - liion.η_self * parameters.Δh) - (p_liion_ch[h] * liion.η_ch + p_liion_dch[h] / liion.η_dch) * parameters.Δh
     # Power balance
     [h in 1:nh], p_net[h] <= p_g_out[h] + p_g_in[h] + p_liion_ch[h] + p_liion_dch[h]
     # Initial and final conditions
     soc_liion[1] == soc_init
     end)

     return m
end

#### Offline functions ####
# Simple
function initialize_controller(ld::Load, pv::Source, liion::Liion, controller::MPCController,
     grid::Grid, ω_optim::Scenarios, parameters::NamedTuple)

     # Parameters
     nh = size(ld.power_E,1) # number of hours
     ny = size(ld.power_E,2) # number of simulation years
     ns = size(ld.power_E,3) # number of scenarios

     # Initialize model
     controller.model = mpc_model(
     ld, pv, liion, controller, grid, parameters)

     # Initialize scenario generation with markov chain
     controller.markovchains = initialize_markovchains(ω_optim)

     # Initialize operation decisions
     controller.u = (
     u_liion =  convert(SharedArray,zeros(nh,ny,ns)),
     )
end

#### Online functions ####
# Simple
function compute_operation_decisions(h::Int64, y::Int64, s::Int64, ld::Load, pv::Source,
     liion::Liion, grid::Grid, controller::MPCController, ω_optim::Scenarios, parameters::NamedTuple)
     # Parameters
     window = h:parameters.Δh:min(parameters.H, h + controller.horizon - 1)
     n_zeros = length(1:parameters.Δh:controller.horizon) - length(window)

     # Compute forecasts and add zeros if needed
     # PV
     pv_fcst = compute_scenario(controller.markovchains.pv_E, pv.power_E[h,y,s] / pv.powerMax[y,s],
      ω_optim.timestamp[h], y, controller.horizon-1)
      # Laod
     ld_E_fcst = compute_scenario(controller.markovchains.ld_E.wk, controller.markovchains.ld_E.wkd,
     ld.power_E[h,y,s], ω_optim.timestamp[h], y, controller.horizon-1)
     # Power net
     power_net_fcst = ld_E_fcst .- pv.powerMax[y,s] .* pv_fcst
     # Grid
     C_grid_in_fcst = vcat(ω_optim.C_grid_in[window,y,s], zeros(n_zeros))
     C_grid_out_fcst = vcat(ω_optim.C_grid_out[window,y,s], zeros(n_zeros))

     # Fix net power
     fix.(controller.model[:p_net], power_net_fcst)

     # Fix initial state
     fix(controller.model[:soc_init], liion.soc[h,y,s] * liion.Erated[y,s])

     # Fix battery capacity
     fix(controller.model[:E_liion], liion.Erated[y,s])

     # Objective function
     @objective(controller.model, Min,
      sum(controller.model[:p_g_in] .* C_grid_in_fcst + controller.model[:p_g_out] .* C_grid_out_fcst) * parameters.Δh)

     # Optimize
     optimize!(controller.model)

     # Operation decision - we only kept the first value
     controller.u.u_liion[h,y,s] = value.(controller.model[:p_liion_ch] .+ controller.model[:p_liion_dch])[1]
end
