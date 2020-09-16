#=
    Model predictive control controller
=#

mutable struct MPCOptions
    solver
    horizon::Int64

    MPCOptions(; solver = CPLEX, horizon = 24) = new(solver, horizon)
end

mutable struct MPC <: AbstractController
    options::MPCOptions
    u::NamedTuple
    model::JuMP.Model
    markovchains
    history::Scenarios
    MPC(; options = MPCOptions()) = new(options)
end

### Model
function build_model!(des::DistributedEnergySystem, controller::MPC, ω_optim::Scenarios)

     # Sets
     nh = controller.options.horizon

     # Model definition
     m = Model(controller.options.solver.Optimizer)
     set_optimizer_attribute(m,"CPX_PARAM_SCRIND", 0)

     # Build model
     if isa(des.liion, Liion)
         # Variables
         @variables(m, begin
         # Operation decisions variables
         p_liion_ch[1:nh] <= 0.
         p_liion_dch[1:nh] >= 0.
         p_g_in[1:nh] >= 0.
         p_g_out[1:nh] <= 0.
         # State variables
         soc_liion[1:nh+1]
         # Variables to be fixed
         p_net[1:nh]
         r_liion
         soc_ini
         end)
         # Constraints
         @constraints(m, begin
         # Power bounds
         [h in 1:nh], p_liion_dch[h] <= des.liion.α_p_dch * r_liion
         [h in 1:nh], p_liion_ch[h] >= -des.liion.α_p_ch * r_liion
         # SoC bounds
         [h in 1:nh+1], soc_liion[h] <= des.liion.α_soc_max * r_liion
         [h in 1:nh+1], soc_liion[h] >= des.liion.α_soc_min * r_liion
         # Dynamics
         [h in 1:nh], soc_liion[h+1] == soc_liion[h] * (1 - des.liion.η_self * des.parameters.Δh) - (p_liion_ch[h] * des.liion.η_ch + p_liion_dch[h] / des.liion.η_dch) * des.parameters.Δh
         # Power balance
         [h in 1:nh], p_net[h] <= p_g_in[h] + p_g_out[h] + p_liion_ch[h] + p_liion_dch[h]
         # Initial and final conditions
         soc_liion[1] == soc_ini
         end)
     end

     # Store model
     controller.model = m
end

### Offline
function initialize_controller!(des::DistributedEnergySystem, controller::MPC, ω_optim::Scenarios)

     # Build model
     build_model!(des, controller, ω_optim)

     # Compute markov chain for scenario generation
     controller.markovchains = compute_markovchains(ω_optim)

     # Store the optimization scenario to the controller history field
     controller.history = ω_optim

     # Preallocate
     preallocate!(controller, des.parameters.nh, des.parameters.ny, des.parameters.ns)

     return controller
end

### Online
function compute_operation_decisions!(h::Int64, y::Int64, s::Int64, des::DistributedEnergySystem, controller::MPC)

     # Parameters
     window = h:min(des.parameters.nh, h + controller.options.horizon - 1)
     n_zeros = controller.options.horizon - length(window)

     # Compute forecasts and add zeros if needed
     pv = compute_scenario(controller.markovchains.pv_E, des.pv.power_E[h,y,s] / des.pv.powerMax[y,s], controller.history.timestamp[h], y, controller.options.horizon - 1)
     ld_E = compute_scenario(controller.markovchains.ld_E.wk, controller.markovchains.ld_E.wkd, des.ld_E.power[h,y,s], controller.history.timestamp[h], y, controller.options.horizon - 1)
     C_grid_in = vcat(controller.history.values.C_grid_in[window,y,s], zeros(n_zeros))
     C_grid_out = vcat(controller.history.values.C_grid_out[window,y,s], zeros(n_zeros))

     # Fix variables
     fix.(controller.model[:p_net], ld_E .- des.pv.powerMax[y,s] .* pv)
     fix(controller.model[:soc_ini], des.liion.soc[h,y,s] * des.liion.Erated[y,s])
     fix(controller.model[:r_liion], des.liion.Erated[y,s])

     # Objective function
     @objective(controller.model, Min, sum(controller.model[:p_g_in] .* C_grid_in + controller.model[:p_g_out] .* C_grid_out) * des.parameters.Δh)

     # Optimize
     optimize!(controller.model)

     # Operation decision - we only keep the first value
     controller.u.liion[h,y,s] = value.(controller.model[:p_liion_ch] .+ controller.model[:p_liion_dch])[1]
end
