#=
    Model predictive control controller
=#

mutable struct OLFCOptions
    solver
    horizon::Int64
    n_state::Int64
    ns::Int64

    OLFCOptions(; solver = CPLEX, horizon = 24, n_state = 5, ns = 10) = new(solver, horizon, n_state, ns)
end

mutable struct OLFC <: AbstractController
    options::OLFCOptions
    u::NamedTuple
    model::JuMP.Model
    markovchains
    history::AbstractScenarios
    OLFC(; options = OLFCOptions()) = new(options)
end

### Model
function build_model(des::DistributedEnergySystem, controller::OLFC)
     # Sets
     nh = controller.options.horizon
     ns = controller.options.ns

     # Initialize
     isa(des.liion, Liion) ? liion = des.liion : liion = Liion()
     isa(des.tes, ThermalSto) ? tes = des.tes : tes = ThermalSto()
     isa(des.h2tank, H2Tank) ? h2tank = des.h2tank : h2tank = H2Tank()
     isa(des.elyz, Electrolyzer) ? elyz = des.elyz : elyz = Electrolyzer()
     isa(des.fc, FuelCell) ? fc = des.fc : fc = FuelCell()
     isa(des.heater, Heater) ? heater = des.heater : heater = Heater()
     isa(des.pv, Source) ? pv = des.pv : pv = Source()
     isa(des.grid, Grid) ? grid = des.grid : grid = Grid()

     # Model definition
     m = Model(controller.options.solver.Optimizer)
     set_optimizer_attribute(m,"CPX_PARAM_SCRIND", 0)

    # Build model

    @variables(m, begin
    # Operation decisions variables
    p_liion_ch[1:nh] <= 0.
    p_liion_dch[1:nh] >= 0.
    p_tes_ch[1:nh] <= 0.
    p_tes_dch[1:nh] >= 0.
    p_h2tank_ch[1:nh] <= 0.
    p_h2tank_dch[1:nh] >= 0.
    p_elyz_E[1:nh] <= 0.
    p_fc_E[1:nh] >= 0.
    p_heater_E[1:nh] <= 0.
    p_g_out[1:nh, 1:ns] <= 0.
    p_g_in[1:nh, 1:ns] >= 0.
    p_net_E[1:nh, 1:ns]
    p_net_H[1:nh, 1:ns]
    # Operation state variables
    soc_liion[1:nh+1]
    soc_tes[1:nh+1]
    soc_h2tank[1:nh+1]
    # Initial states to be fixed
    soc_liion_ini
    soc_tes_ini
    soc_h2tank_ini
    # Investment variables to be fixed
    r_liion == 0
    r_tes == 0
    r_h2tank == 0
    r_elyz == 0
    r_fc == 0
    end)
    @constraints(m, begin
    # Power bounds
    [h in 1:nh], p_liion_dch[h] <= liion.α_p_dch * r_liion
    [h in 1:nh], p_liion_ch[h] >= -liion.α_p_ch * r_liion
    [h in 1:nh], p_tes_dch[h] <= tes.α_p_dch * r_tes
    [h in 1:nh], p_tes_ch[h] >= -tes.α_p_ch * r_tes
    [h in 1:nh], p_h2tank_dch[h] <= h2tank.α_p_dch * r_h2tank
    [h in 1:nh], p_h2tank_ch[h] >= -h2tank.α_p_ch * r_h2tank
    [h in 1:nh], p_elyz_E[h] >= -r_elyz
    [h in 1:nh], p_fc_E[h] <= r_fc
    [h in 1:nh], p_heater_E[h] >= -heater.powerMax_ini
    [h in 1:nh, s in 1:ns], p_g_in[h,s] <= grid.powerMax
    # SoC bounds
    [h in 1:nh+1], soc_liion[h] <= liion.α_soc_max * r_liion
    [h in 1:nh+1], soc_liion[h] >= liion.α_soc_min * r_liion
    [h in 1:nh+1], soc_tes[h] <= tes.α_soc_max * r_tes
    [h in 1:nh+1], soc_tes[h] >= tes.α_soc_min * r_tes
    [h in 1:nh+1], soc_h2tank[h] <= h2tank.α_soc_max * r_h2tank
    [h in 1:nh+1], soc_h2tank[h] >= h2tank.α_soc_min * r_h2tank
    # State dynamics
    [h in 1:nh], soc_liion[h+1] == soc_liion[h] * (1. - liion.η_self * des.parameters.Δh) - (p_liion_ch[h] * liion.η_ch + p_liion_dch[h] / liion.η_dch) * des.parameters.Δh
    [h in 1:nh], soc_tes[h+1] == soc_tes[h] * (1. - tes.η_self * des.parameters.Δh) - (p_tes_ch[h] * tes.η_ch + p_tes_dch[h] / tes.η_dch) * des.parameters.Δh
    [h in 1:nh], soc_h2tank[h+1] == soc_h2tank[h] * (1. - h2tank.η_self * des.parameters.Δh) - (p_h2tank_ch[h] * h2tank.η_ch + p_h2tank_dch[h] / h2tank.η_dch) * des.parameters.Δh
    # Initial and final conditions
    soc_liion[1] == soc_liion_ini
    soc_tes[1] == soc_tes_ini
    soc_h2tank[1] == soc_h2tank_ini
    # Power balances
    [h in 1:nh, s in 1:ns], p_net_E[h,s] <= p_liion_ch[h] + p_liion_dch[h] + p_elyz_E[h] + p_fc_E[h] + p_heater_E[h] + p_g_in[h,s] + p_g_out[h,s]
    [h in 1:nh, s in 1:ns], p_net_H[h,s] <= p_tes_ch[h] + p_tes_dch[h] - elyz.η_E_H .* p_elyz_E[h] + fc.η_H2_H / fc.η_H2_E .* p_fc_E[h] - heater.η_E_H .* p_heater_E[h]
    [h in 1:nh, s in 1:ns], 0. == p_h2tank_ch[h] + p_h2tank_dch[h] - elyz.η_E_H2 * p_elyz_E[h] - p_fc_E[h] / fc.η_H2_E
    end)

    return m
end
### Offline
function initialize_controller!(des::DistributedEnergySystem, controller::OLFC, ω::AbstractScenarios)

     # Build model
     controller.model = build_model(des, controller)

     # Scenario reduciton
     ω_markov = reduce(ω, 1:des.parameters.nh, 1:des.parameters.ny, 1)

     # Compute markov chain for scenario generation TODO pas necessaire acr toujours ld_E dans scenarios
     if isa(des.ld_H, Load)
         mc_wkd, mc_wk = compute_markovchains(ω_markov.pv, ω_markov.ld_E, ω_markov.ld_H, n_state = controller.options.n_state)
     else
         mc_wkd, mc_wk = compute_markovchains(ω_markov.pv, ω_markov.ld_E, n_state = controller.options.n_state)
     end
     controller.markovchains = (wk = mc_wk, wkd = mc_wkd)

     # Preallocate
     preallocate!(controller, des.parameters.nh, des.parameters.ny, des.parameters.ns)

     return controller
end

### Online
function compute_operation_decisions!(h::Int64, y::Int64, s::Int64, des::DistributedEnergySystem, controller::OLFC)
    # Parameters
    window = h:min(des.parameters.nh, h + controller.options.horizon - 1)
    n_zeros = controller.options.horizon - length(window)
    isa(des.ld_H, Load) ? s0 = [des.pv.power_E[h,y,s] / des.pv.powerMax[y,s], des.ld_E.power[h,y,s], des.ld_H.power[h,y,s]] : s0 = [des.pv.power_E[h,y,s] / des.pv.powerMax[y,s], des.ld_E.power[h,y,s]]
    t0 = des.pv.timestamp[h,y,s]

    # Compute forecast
    forecast = Genesys.compute_scenarios(controller.markovchains.wk, controller.markovchains.wkd, s0, t0, controller.options.horizon, ny = controller.options.ns)
    cost_in = vcat(des.grid.cost_in[window,y,s], zeros(n_zeros)) ; cost_out = vcat(des.grid.cost_out[window,y,s], zeros(n_zeros))

    # Fix forecast and state variables
    fix_variables!(h, y, s, des, controller, forecast)

    # Objective function TODO ajouter correctement les probas !
    @objective(controller.model, Min, sum(controller.model[:p_g_in] .* cost_in .+ controller.model[:p_g_out] .* cost_out) * des.parameters.Δh / controller.options.ns)

    # Optimize
    optimize!(controller.model)

    # Operation decision - we only keep the first value
    controller.u.liion[h,y,s] = value.(controller.model[:p_liion_ch] + controller.model[:p_liion_dch])[1]
    controller.u.tes[h,y,s] = value.(controller.model[:p_tes_ch] + controller.model[:p_tes_dch])[1]
    controller.u.h2tank[h,y,s] = value.(controller.model[:p_h2tank_ch] + controller.model[:p_h2tank_dch])[1]
    controller.u.elyz[h,y,s] = value.(controller.model[:p_elyz_E])[1]
    controller.u.fc[h,y,s] = value.(controller.model[:p_fc_E])[1]
    controller.u.heater[h,y,s] = value.(controller.model[:p_heater_E])[1]
end

### Utils
function fix_variables!(h::Int64, y::Int64, s::Int64, des::DistributedEnergySystem, controller::OLFC, forecast)
    # Fix forecast and state variables
    isa(des.ld_E, Load) ? fix.(controller.model[:p_net_E], forecast[2].power .- des.pv.powerMax[y,s] .* forecast[1].power) : fix.(controller.model[:p_net_E], 0.)
    isa(des.ld_H, Load) ? fix.(controller.model[:p_net_H], forecast[3].power) : fix.(controller.model[:p_net_H], 0.)

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
end
