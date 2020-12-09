# Anticipative controller

mutable struct AnticipativeOptions
    solver

    AnticipativeOptions(; solver = CPLEX) = new(solver)
end

mutable struct Anticipative <: AbstractController
    options::AnticipativeOptions
    pv::Float64
    liion::Float64
    tes::Float64
    h2tank::Float64
    elyz::Float64
    fc::Float64
    u::NamedTuple

    Anticipative(; options = AnticipativeOptions(),
                   pv = 0.,
                   liion = 0.,
                   tes = 0.,
                   h2tank = 0.,
                   elyz = 0.,
                   fc = 0.) =
                   new(options, pv, liion, tes, h2tank, elyz, fc)
end

### Models
function build_model(des::DistributedEnergySystem, controller::Anticipative, ω::Scenarios)
    # Sets
    nh = size(ω.ld_E.power,1) # Number of hours

    # Initialize
    isa(des.ld_E, Load) ? ld_E = ω.ld_E.power : ld_E = zeros(nh)
    isa(des.ld_H, Load) ? ld_H = ω.ld_H.power : ld_H = zeros(nh)
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
    # Variables
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
    p_g_out[1:nh] <= 0.
    p_g_in[1:nh] >= 0.
    # Operation state variables
    soc_liion[1:nh+1]
    soc_tes[1:nh+1]
    soc_h2tank[1:nh+1]
    end)
    # Constraints
    @constraints(m, begin
    # Power bounds
    [h in 1:nh], p_liion_dch[h] <= liion.α_p_dch * controller.liion
    [h in 1:nh], p_liion_ch[h] >= -liion.α_p_ch * controller.liion
    [h in 1:nh], p_tes_dch[h] <= tes.α_p_dch * controller.tes
    [h in 1:nh], p_tes_ch[h] >= -tes.α_p_ch * controller.tes
    [h in 1:nh], p_h2tank_dch[h] <= h2tank.α_p_dch * controller.h2tank
    [h in 1:nh], p_h2tank_ch[h] >= -h2tank.α_p_ch * controller.h2tank
    [h in 1:nh], p_elyz_E[h] >= -controller.elyz
    [h in 1:nh], p_fc_E[h] <= controller.fc
    [h in 1:nh], p_heater_E[h] >= -heater.powerMax_ini
    [h in 1:nh], p_g_in[h] <= grid.powerMax
    [h in 1:nh], p_g_out[h] >= -grid.powerMax
    # SoC bounds
    [h in 1:nh+1], soc_liion[h] <= liion.α_soc_max * controller.liion
    [h in 1:nh+1], soc_liion[h] >= liion.α_soc_min * controller.liion
    [h in 1:nh+1], soc_tes[h] <= tes.α_soc_max * controller.tes
    [h in 1:nh+1], soc_tes[h] >= tes.α_soc_min * controller.tes
    [h in 1:nh+1], soc_h2tank[h] <= h2tank.α_soc_max * controller.h2tank
    [h in 1:nh+1], soc_h2tank[h] >= h2tank.α_soc_min * controller.h2tank
    # State dynamics
    [h in 1:nh], soc_liion[h+1] == soc_liion[h] * (1. - liion.η_self * des.parameters.Δh) - (p_liion_ch[h] * liion.η_ch + p_liion_dch[h] / liion.η_dch) * des.parameters.Δh
    [h in 1:nh], soc_tes[h+1] == soc_tes[h] * (1. - tes.η_self * des.parameters.Δh) - (p_tes_ch[h] * tes.η_ch + p_tes_dch[h] / tes.η_dch) * des.parameters.Δh
    [h in 1:nh], soc_h2tank[h+1] == soc_h2tank[h] * (1. - h2tank.η_self * des.parameters.Δh) - (p_h2tank_ch[h] * h2tank.η_ch + p_h2tank_dch[h] / h2tank.η_dch) * des.parameters.Δh
    # Initial and final states
    soc_liion[1] == liion.soc_ini * controller.liion
    soc_liion[end] >= soc_liion[1]
    soc_tes[1] == tes.soc_ini * controller.tes
    soc_tes[end] >= soc_tes[1]
    soc_h2tank[1] == h2tank.soc_ini * controller.h2tank
    soc_h2tank[end] >= soc_h2tank[1]
    # Power balances
    [h in 1:nh], ld_E[h] <= controller.pv * ω.pv.power[h] + p_liion_ch[h] + p_liion_dch[h] + p_elyz_E[h] + p_fc_E[h] + p_heater_E[h] + p_g_in[h] + p_g_out[h]
    [h in 1:nh], ld_H[h] <= p_tes_ch[h]  + p_tes_dch[h] - elyz.η_E_H * p_elyz_E[h] + fc.η_H2_H / fc.η_H2_E * p_fc_E[h] - heater.η_E_H * p_heater_E[h]
    [h in 1:nh], 0. == p_h2tank_ch[h] + p_h2tank_dch[h] - elyz.η_E_H2 * p_elyz_E[h] - p_fc_E[h] / fc.η_H2_E
    end)

    # Share of renewables constraint - max(0, share)
    @expression(m, share_constraint, sum(p_g_in[h] - (1. - des.parameters.τ_share) * (ld_E[h] + ld_H[h] / heater.η_E_H) for h in 1:nh))
    @variable(m, aux >= 0.)
    @constraint(m, aux >= share_constraint)

    # OPEX
    @expression(m, opex, sum((p_g_in[h] * ω.grid.cost_in[h] + p_g_out[h] * ω.grid.cost_out[h]) * des.parameters.Δh  for h in 1:nh))

    # Objective
    λ = 1e6
    @objective(m, Min, opex + λ * aux)

    return m
end

### Offline
function initialize_controller!(des::DistributedEnergySystem, controller::Anticipative, ω::Scenarios)
    # Preallocate
    preallocate!(controller, des.parameters.nh, des.parameters.ny, des.parameters.ns)

    for y in 2:des.parameters.ny, s in 1:des.parameters.ns
        # Scenario reduction
        ω_reduced, _ = reduce(ManualReducer(h = 1:des.parameters.nh, y = y, s = s), ω)
        # Build model
        model = build_model(des, controller, ω_reduced)
        # Optimize
        optimize!(model)
        # Assign controller values
        controller.u.liion[:,y,s] .= value.(model[:p_liion_ch] + model[:p_liion_dch])
        controller.u.tes[:,y,s] .= value.(model[:p_tes_ch] + model[:p_tes_dch])
        controller.u.h2tank[:,y,s] .= value.(model[:p_h2tank_ch] + model[:p_h2tank_dch])
        controller.u.elyz[:,y,s] .= value.(model[:p_elyz_E])
        controller.u.fc[:,y,s] .= value.(model[:p_fc_E])
        controller.u.heater[:,y,s] .= value.(model[:p_heater_E])
     end

     return controller
end

### Online
function compute_operation_decisions!(h::Int64, y::Int64, s::Int64, des::DistributedEnergySystem, controller::Anticipative)
    return controller
end
