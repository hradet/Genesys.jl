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
    history::AbstractScenarios
    MPC(; options = MPCOptions()) = new(options)
end

### Model
function build_model(des::DistributedEnergySystem, controller::MPC, ω::AbstractScenarios)

     # Sets
     nh = controller.options.horizon

     # Model definition
     m = Model(controller.options.solver.Optimizer)
     set_optimizer_attribute(m,"CPX_PARAM_SCRIND", 0)

    # Build model
    ### Initilize power balances
    if isa(des.ld_E, Load)
        # Variable to be fixed
        @variables(m, begin
        p_net_E[1:nh]
        p_ld_E[1:nh] >= 0.
        end)
        @constraints(m, begin
        [h in 1:nh], p_ld_E[h] >= p_net_E[h]
        end)
        power_balance_E = @expression(m, - p_net_E)
    end
    if isa(des.ld_H, Load)
        # Variable to be fixed
        @variables(m, begin
        p_net_H[1:nh]
        p_ld_H[1:nh] >= 0.
        end)
        @constraints(m, begin
        [h in 1:nh], p_ld_H[h] >= p_net_H[h]
        end)
        power_balance_H = @expression(m, - p_net_H)
    end
    power_balance_H2 = AffExpr.(zeros(nh))

    ### Storages
    if isa(des.liion, Liion)
        # Variables
        @variables(m, begin
        # Operation decisions variables
        p_liion_ch[1:nh] <= 0.
        p_liion_dch[1:nh] >= 0.
        # Operation state variables
        soc_liion[1:nh+1]
        # To be fixed
        r_liion
        soc_liion_ini
        end)
        # Constraints
        @constraints(m, begin
        # Power bounds
        [h in 1:nh], p_liion_dch[h] <= des.liion.α_p_dch * r_liion
        [h in 1:nh], p_liion_ch[h] >= -des.liion.α_p_ch * r_liion
        # SoC bounds
        [h in 1:nh+1], soc_liion[h] <= des.liion.α_soc_max * r_liion
        [h in 1:nh+1], soc_liion[h] >= des.liion.α_soc_min * r_liion
        # State dynamics
        [h in 1:nh], soc_liion[h+1] == soc_liion[h] * (1. - des.liion.η_self * des.parameters.Δh) - (p_liion_ch[h] * des.liion.η_ch + p_liion_dch[h] / des.liion.η_dch) * des.parameters.Δh
        # Initial and final conditions
        soc_liion[1] == soc_liion_ini
        soc_liion[nh] >= soc_liion_ini
        end)
        # Power balance
        isa(des.ld_E, Load) ? add_to_expression!.(power_balance_E, p_liion_ch .+ p_liion_dch) : nothing
    end
    if isa(des.tes, ThermalSto)
        # Variables
        @variables(m, begin
        # Operation decisions variables
        p_tes_ch[1:nh] <= 0.
        p_tes_dch[1:nh] >= 0.
        # Operation state variables
        soc_tes[1:nh+1]
        # To be fixed
        r_tes
        soc_tes_ini
        end)
        # Constraints
        @constraints(m, begin
        # Power bounds
        [h in 1:nh], p_tes_dch[h] <= des.tes.α_p_dch * r_tes
        [h in 1:nh], p_tes_ch[h] >= -des.tes.α_p_ch * r_tes
        # SoC bounds
        [h in 1:nh+1], soc_tes[h] <= des.tes.α_soc_max * r_tes
        [h in 1:nh+1], soc_tes[h] >= des.tes.α_soc_min * r_tes
        # State dynamics
        [h in 1:nh], soc_tes[h+1] == soc_tes[h] * (1. - des.tes.η_self * des.parameters.Δh) - (p_tes_ch[h] * des.tes.η_ch + p_tes_dch[h] / des.tes.η_dch) * des.parameters.Δh
        # Initial and final conditions
        soc_tes[1] == soc_tes_ini
        soc_tes[nh] >= soc_tes_ini
        end)
        # Power balance
        isa(des.ld_H, Load) ? add_to_expression!.(power_balance_H, p_tes_ch .+ p_tes_dch) : nothing
     end
     if isa(des.h2tank, H2Tank)
        # Variables
        @variables(m, begin
        # Operation decisions variables
        p_h2tank_ch[1:nh] <= 0.
        p_h2tank_dch[1:nh] >= 0.
        # Operation state variables
        soc_h2tank[1:nh+1]
        # To be fixed
        r_h2tank
        soc_h2tank_ini
        end)
        # Constraints
        @constraints(m, begin
        # Power bounds
        [h in 1:nh], p_h2tank_dch[h] <= des.h2tank.α_p_dch * r_h2tank
        [h in 1:nh], p_h2tank_ch[h] >= -des.h2tank.α_p_ch * r_h2tank
        # SoC bounds
        [h in 1:nh+1], soc_h2tank[h] <= des.h2tank.α_soc_max * r_h2tank
        [h in 1:nh+1], soc_h2tank[h] >= des.h2tank.α_soc_min * r_h2tank
        # State dynamics
        [h in 1:nh], soc_h2tank[h+1] == soc_h2tank[h] * (1. - des.h2tank.η_self * des.parameters.Δh) - (p_h2tank_ch[h] * des.h2tank.η_ch + p_h2tank_dch[h] / des.h2tank.η_dch) * des.parameters.Δh
        # Initial and final conditions
        soc_h2tank[1] == soc_h2tank_ini
        soc_h2tank[nh] >= soc_h2tank_ini
        end)
        # Power balances
        add_to_expression!.(power_balance_H2, p_h2tank_ch .+ p_h2tank_dch)
    end
    ### Converters
    if isa(des.elyz, Electrolyzer)
        # Variables
        @variables(m, begin
        # Operation decisions variables
        p_elyz_E[1:nh] <= 0.
        # To be fixed
        r_elyz
        end)
        # Constraints
        @constraints(m, begin
        # Power bounds
        [h in 1:nh], p_elyz_E[h] >= -r_elyz
        end)
        # Power balance
        isa(des.ld_E, Load) ? add_to_expression!.(power_balance_E, p_elyz_E) : nothing
        isa(des.ld_H, Load) ? add_to_expression!.(power_balance_H, - des.elyz.η_E_H .* p_elyz_E) : nothing
        add_to_expression!.(power_balance_H2, - des.elyz.η_E_H2 * p_elyz_E)
    end
    if isa(des.fc, FuelCell)
        # Variables
        @variables(m, begin
        # Operation decisions variables
        p_fc_E[1:nh] >= 0.
        # To be fixed
        r_fc
        end)
        # Constraints
        @constraints(m, begin
        # Power bounds
        [h in 1:nh], p_fc_E[h] <= r_fc
        end)
        # Power balance
        isa(des.ld_E, Load) ? add_to_expression!.(power_balance_E, p_fc_E) : nothing
        isa(des.ld_H, Load) ? add_to_expression!.(power_balance_H, des.fc.η_H2_H / des.fc.η_H2_E .* p_fc_E) : nothing
        add_to_expression!.(power_balance_H2, - p_fc_E / des.fc.η_H2_E)
    end
    if isa(des.heater, Heater)
        # Variables
        @variables(m, begin
        # Operation decisions variables
        p_heater_E[1:nh] <= 0.
        end)
        # Constraints
        @constraints(m, begin
        # Power bounds
        [h in 1:nh], p_heater_E[h] >= -des.heater.powerMax[1]
        end)
        # Power balance
        isa(des.ld_E, Load) ? add_to_expression!.(power_balance_E, p_heater_E) : nothing
        isa(des.ld_H, Load) ? add_to_expression!.(power_balance_H, - des.heater.η_E_H .* p_heater_E) : nothing
    end
    ### Grid
    if isa(des.grid, Grid)
        # Variables
        @variables(m, begin
        # Operation decisions variables
        p_g_out[1:nh] <= 0.
        p_g_in[1:nh] >= 0.
        end)
        # Constraints
        @constraints(m, begin
        # Power bounds
        [h in 1:nh], p_g_in[h] <= des.grid.powerMax
        end)
        # Power balance
        isa(des.ld_E, Load) ? add_to_expression!.(power_balance_E, p_g_in .+ p_g_out) : nothing
    end
    ### Power balance constraints
    isa(des.ld_E, Load) ? @constraint(m, power_balance_E .>= 0.) : nothing
    isa(des.ld_H, Load) ? @constraint(m, power_balance_H .>= 0.) : nothing
    @constraint(m, power_balance_H2 .== 0.)

    return m
end

### Offline
function initialize_controller!(des::DistributedEnergySystem, controller::MPC, ω::AbstractScenarios)

     # Build model
     controller.model = build_model(des, controller, ω)

     # Compute markov chain for scenario generation TODO gérer la chaleur !
     controller.markovchains = compute_markovchains(ω.pv, ω.ld_E, ω.ld_H)

     # Store the optimization scenario to the controller history field
     controller.history = ω

     # Preallocate
     preallocate!(controller, des.parameters.nh, des.parameters.ny, des.parameters.ns)

     return controller
end

### Online
function compute_operation_decisions!(h::Int64, y::Int64, s::Int64, des::DistributedEnergySystem, controller::MPC)

     # Parameters
     window = h:min(des.parameters.nh, h + controller.options.horizon - 1)
     n_zeros = controller.options.horizon - length(window)
     s0 = [des.pv.power_E[h,y,s], des.ld_E.power_E[h,y,s], des.ld_E.power_H[h,y,s]]
     t0 = des.pv.timestamp[h,y,s]

     # Compute forecasts
     if isa(des.ld_E, Load)
         pv, ld_E, ld_H = compute_scenario(controller.markovchains, s0, t0, controller.options.horizon - 1)
         fix.(controller.model[:p_net_E], ld_E.power .- des.pv.powerMax[y,s] .* pv.power)
         fix.(controller.model[:p_net_H], ld_H)
     end
     # Fix state variables
     if isa(des.liion, Liion)
         fix(controller.model[:soc_liion_ini], des.liion.soc[h,y,s] * des.liion.Erated[y,s])
         fix(controller.model[:r_liion], des.liion.Erated[y,s])
     end
     if isa(des.tes, ThermalSto)
         fix(controller.model[:soc_tes_ini], des.tes.soc[h,y,s] * des.tes.Erated[y,s])
         fix(controller.model[:r_tes], des.tes.Erated[y,s])
     end
     if isa(des.h2tank, H2Tank)
         fix(controller.model[:soc_h2tank_ini], des.h2tank.soc[h,y,s] * des.h2tank.Erated[y,s])
         fix(controller.model[:r_h2tank], des.h2tank.Erated[y,s])
     end
     isa(des.elyz, Electrolyzer) ? fix(controller.model[:r_elyz], des.elyz.powerMax[y,s]) : nothing
     isa(des.fc, FuelCell) ? fix(controller.model[:r_fc], des.fc.powerMax[y,s]) : nothing
     cost_in = vcat(controller.history.grid.cost_in[window,y,s], zeros(n_zeros))
     cost_out = vcat(controller.history.grid.cost_out[window,y,s], zeros(n_zeros))

     # Objective function
     @objective(controller.model, Min, sum(controller.model[:p_g_in] .* cost_in + controller.model[:p_g_out] .* cost_out) * des.parameters.Δh)

     # Optimize
     optimize!(controller.model)

     # Operation decision - we only keep the first value
     isa(des.liion, Liion) ? controller.u.liion[h,y,s] = value.(controller.model[:p_liion_ch] .+ controller.model[:p_liion_dch])[1] : nothing
     isa(des.tes, ThermalSto) ? controller.u.tes[h,y,s] = value.(controller.model[:p_tes_ch] .+ controller.model[:p_tes_dch])[1] : nothing
     isa(des.h2tank, H2Tank) ? controller.u.h2tank[h,y,s] = value.(controller.model[:p_h2tank_ch] .+ controller.model[:p_h2tank_dch])[1] : nothing
     isa(des.elyz, Electrolyzer) ? controller.u.elyz[h,y,s] = value.(controller.model[:p_elyz_E])[1] : nothing
     isa(des.fc, FuelCell) ? controller.u.fc[h,y,s] = value.(controller.model[:p_fc_E])[1] : nothing
     isa(des.heater, Heater) ? controller.u.heater[h,y,s] = value.(controller.model[:p_heater_E])[1] : nothing

end
