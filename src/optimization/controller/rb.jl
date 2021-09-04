#=
    Rule based controller
=#
mutable struct RBCOptions
    policy_selection::Int64

    RBCOptions(; policy_selection = 1) = new(policy_selection)
end

mutable struct RBC <: AbstractController
    options::RBCOptions
    decisions::NamedTuple
    history::AbstractScenarios

    RBC(; options = RBCOptions()) = new(options)
end

### Policies
function π_1(h::Int64, y::Int64, s::Int64, mg::Microgrid, controller::RBC)
    # Utils to simplify the writting
    Δh = mg.parameters.Δh
    liion, tes, h2tank = mg.storages[1], mg.storages[2], mg.storages[3]
    heater, elyz, fc = mg.converters[1], mg.converters[2], mg.converters[3]

    # Net power elec
    p_net_E = mg.demands[1].carrier.out[h,y,s] - mg.generations[1].carrier.in[h,y,s]

    if p_net_E < 0
        # Liion
        power_in = max(min(p_net_E, liion.α_p_dch * liion.Erated[y,s], liion.soh[h,y,s] * liion.Erated[y,s] / Δh, liion.η_dch * (liion.soc[h,y,s] * (1. - liion.η_self * Δh) - liion.α_soc_min) * liion.Erated[y,s] / Δh), 0.)
        power_out = min(max(p_net_E, -liion.α_p_ch * liion.Erated[y,s], -liion.soh[h,y,s] * liion.Erated[y,s] / Δh, (liion.soc[h,y,s] * (1. - liion.η_self * Δh) - liion.α_soc_max) * liion.Erated[y,s] / Δh / liion.η_ch), 0.)
        u_liion = power_in + power_out
        # Elyz
        elyz.α_p * elyz.powerMax[y,s] >= p_net_E - u_liion && elyz.soh[h,y,s] * elyz.nHoursMax / Δh > 1. ? u_elyz_E = max(p_net_E - u_liion, -elyz.powerMax[y,s]) : u_elyz_E = 0.
        elyz_H2 = - u_elyz_E * elyz.η_E_H2
        elyz_H = - u_elyz_E * elyz.η_E_H
        # H2 tank
        power_in = max(min(- elyz_H2, h2tank.α_p_dch * h2tank.Erated[y,s], h2tank.η_dch * (h2tank.soc[h,y,s] * (1. - h2tank.η_self * Δh) - h2tank.α_soc_min) * h2tank.Erated[y,s] / Δh), 0.)
        power_out = min(max(- elyz_H2, -h2tank.α_p_ch * h2tank.Erated[y,s], (h2tank.soc[h,y,s] * (1. - h2tank.η_self * Δh) - h2tank.α_soc_max) * h2tank.Erated[y,s] / Δh / h2tank.η_ch), 0.)
        u_h2tank = power_in + power_out
        # Test H2
        elyz_H2 == - u_h2tank ? nothing : u_elyz_E = elyz_H = elyz_H2 = u_h2tank = 0.
        # FC
        u_fc_E, fc_H, fc_H2 = 0., 0., 0.
        # Heater
        u_heater_E = min(max(p_net_E - u_liion - u_elyz_E, -heater.powerMax[y,s]), 0.)
        heater_H = - heater.η_E_H * u_heater_E
    else
        # Liion
        power_in = max(min(p_net_E, liion.α_p_dch * liion.Erated[y,s], liion.soh[h,y,s] * liion.Erated[y,s] / Δh, liion.η_dch * (liion.soc[h,y,s] * (1. - liion.η_self * Δh) - liion.α_soc_min) * liion.Erated[y,s] / Δh), 0.)
        power_out = min(max(p_net_E, -liion.α_p_ch * liion.Erated[y,s], -liion.soh[h,y,s] * liion.Erated[y,s] / Δh, (liion.soc[h,y,s] * (1. - liion.η_self * Δh) - liion.α_soc_max) * liion.Erated[y,s] / Δh / liion.η_ch), 0.)
        u_liion = power_in + power_out
        # FC
        fc.α_p * fc.powerMax[y,s] <= p_net_E - u_liion && fc.soh[h,y,s] * fc.nHoursMax / Δh > 1. ? u_fc_E = min(p_net_E - u_liion, fc.powerMax[y,s]) : u_fc_E = 0.
        # Power conversion
        fc_H2 = - u_fc_E / fc.η_H2_E
        fc_H = - fc_H2 * fc.η_H2_H
        # H2 tank
        power_in = max(min(- fc_H2, h2tank.α_p_dch * h2tank.Erated[y,s], h2tank.η_dch * (h2tank.soc[h,y,s] * (1. - h2tank.η_self * Δh) - h2tank.α_soc_min) * h2tank.Erated[y,s] / Δh), 0.)
        power_out = min(max(- fc_H2, -h2tank.α_p_ch * h2tank.Erated[y,s], (h2tank.soc[h,y,s] * (1. - h2tank.η_self * Δh) - h2tank.α_soc_max) * h2tank.Erated[y,s] / Δh / h2tank.η_ch), 0.)
        u_h2tank = power_in + power_out
        # Test H2
        fc_H2 == - u_h2tank ? nothing : u_fc_E = fc_H = fc_H2 = u_h2tank = 0.
        # Elyz
        u_elyz_E, elyz_H, elyz_H2 = 0., 0., 0.
        # Heater
        u_heater_E, heater_H = 0., 0.
    end

    # Net heating power post H2
    p_net_H = mg.demands[2].carrier.out[h,y,s] - fc_H - elyz_H - heater_H

    if p_net_H < 0
        # TES
        power_in = max(min(p_net_H, tes.α_p_dch * tes.Erated[y,s], tes.η_dch * (tes.soc[h,y,s] * (1. - tes.η_self * Δh) - tes.α_soc_min) * tes.Erated[y,s] / Δh), 0.)
        power_out = min(max(p_net_H, -tes.α_p_ch * tes.Erated[y,s], (tes.soc[h,y,s] * (1. - tes.η_self * Δh) - tes.α_soc_max) * tes.Erated[y,s] / Δh / tes.η_ch), 0.)
        u_tes = power_in + power_out
        # Heater
        _u_heater_E = 0.
    else
        # TES
        power_in = max(min(p_net_H, tes.α_p_dch * tes.Erated[y,s], tes.η_dch * (tes.soc[h,y,s] * (1. - tes.η_self * Δh) - tes.α_soc_min) * tes.Erated[y,s] / Δh), 0.)
        power_out = min(max(p_net_H, -tes.α_p_ch * tes.Erated[y,s], (tes.soc[h,y,s] * (1. - tes.η_self * Δh) - tes.α_soc_max) * tes.Erated[y,s] / Δh / tes.η_ch), 0.)
        u_tes = power_in + power_out
        # Heater
        _u_heater_E = min(max(- (p_net_H - u_tes), -heater.powerMax[y,s]), 0.)
    end

    # Store values
    controller.decisions.storages[1][h,y,s] = u_liion
    controller.decisions.storages[2][h,y,s] = u_tes
    controller.decisions.storages[3][h,y,s] = u_h2tank
    controller.decisions.converters[1][h,y,s]  = u_heater_E + _u_heater_E
    controller.decisions.converters[2][h,y,s] = u_elyz_E
    controller.decisions.converters[3][h,y,s] = u_fc_E
end
function π_2(h::Int64, y::Int64, s::Int64, mg::Microgrid, controller::RBC)
    controller.decisions.storages[1][h,y,s] = mg.demands[1].carrier.out[h,y,s] - mg.generations[1].carrier.in[h,y,s]
end
function π_3(h::Int64, y::Int64, s::Int64, mg::Microgrid, controller::RBC)
    # Compute the heater electrical power based on the simple model
    controller.decisions.converters[1][h,y,s] = - max(min(mg.demands[2].carrier.out[h,y,s] / mg.converters[1].η_E_H, mg.converters[1].powerMax[y,s]), 0.)
    # Compute the liion decision from the power balance
    controller.decisions.storages[1][h,y,s] = mg.demands[1].carrier.out[h,y,s] - mg.generations[1].carrier.in[h,y,s] - controller.decisions.converters[1][h,y,s]
end


### Offline
function initialize_controller!(mg::Microgrid, controller::RBC, ω::AbstractScenarios)
    # Preallocation
    preallocate!(mg, controller)

    return controller
end

### Online
function compute_operation_decisions!(h::Int64, y::Int64, s::Int64, mg::Microgrid, controller::RBC)
    # Chose policy
    if controller.options.policy_selection == 1
        return π_1(h, y, s, mg, controller)
    elseif controller.options.policy_selection == 2
        return π_2(h, y, s, mg, controller)
    elseif controller.options.policy_selection == 3
        return π_3(h, y, s, mg, controller)
    else
        println("Policy not defined !")
    end
end
