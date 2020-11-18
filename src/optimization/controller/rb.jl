#=
    Rule based controller
=#
mutable struct RBCOptions
    β_min_tes::Float64
    β_max_tes::Float64
    β_min_fc::Float64
    β_max_fc::Float64
    β_min_elyz::Float64
    β_max_elyz::Float64
    optim_threshold::Float64

    RBCOptions(; β_min_tes = 0.2,
                 β_max_tes = 0.9,
                 β_min_fc = 0.25,
                 β_max_fc = 0.3,
                 β_min_elyz = 0.45,
                 β_max_elyz = 0.5,
                 optim_threshold = false) =
                 new(β_min_tes, β_max_tes, β_min_fc, β_max_fc, β_min_elyz, β_max_elyz,optim_threshold)
end

mutable struct RBC <: AbstractController
    options::RBCOptions
    u::NamedTuple
    history::AbstractScenarios

    RBC(; options = RBCOptions()) = new(options)
end

### Policies
function π_1(h::Int64, y::Int64, s::Int64, des::DistributedEnergySystem, controller::RBC)
    # Control parameters
    ϵ = 1e-2

    # Net power elec
    p_net_E = des.ld_E.power[h,y,s] - des.pv.power_E[h,y,s]

    # H2 - Ulleberg, 2004
    if p_net_E <= 0
        # FC is off
        u_fc_E, p_fc_H = 0., 0.
        # Strategy for Elyz
        if des.liion.soc[h,y,s] > controller.options.β_max_elyz
            u_elyz_E = min(max(p_net_E, -des.elyz.powerMax[y,s]), -des.elyz.α_p * des.elyz.powerMax[y,s])
            p_elyz_H = - u_elyz_E * des.elyz.η_E_H
        elseif des.liion.soc[h,y,s] > controller.options.β_min_elyz && h !=1 && !isapprox(des.elyz.power_E[h-1,y,s], 0., atol = ϵ)
            u_elyz_E = min(max(p_net_E, -des.elyz.powerMax[y,s]), -des.elyz.α_p * des.elyz.powerMax[y,s])
            p_elyz_H = - u_elyz_E * des.elyz.η_E_H
        else
            u_elyz_E, p_elyz_H = 0., 0.
        end
    else
        # Elyz is off
        u_elyz_E, p_elyz_H = 0., 0.
        # Strategy for FC
        if des.liion.soc[h,y,s] < controller.options.β_min_fc
            u_fc_E = min(max(p_net_E, des.fc.α_p * des.fc.powerMax[y,s]), des.fc.powerMax[y,s])
            p_fc_H = u_fc_E * des.fc.η_H2_H / des.fc.η_H2_E
        elseif des.liion.soc[h,y,s] < controller.options.β_max_fc && h !=1 && !isapprox(des.fc.power_E[h-1,y,s], 0., atol = ϵ)
            u_fc_E = min(max(p_net_E, des.fc.α_p * des.fc.powerMax[y,s]), des.fc.powerMax[y,s])
            p_fc_H = u_fc_E * des.fc.η_H2_H / des.fc.η_H2_E
        else
            u_fc_E, p_fc_H = 0., 0.
        end
    end

    # Net power elec post H2
    p_net_E += - u_elyz_E - u_fc_E

    # Net power heating post H2
    p_net_H = des.ld_H.power[h,y,s] - p_fc_H - p_elyz_H

    # Heater
    if p_net_H >= 0.
        if des.tes.soc[h,y,s] < controller.options.β_min_tes
            p_heater_H = p_net_H
            u_heater_E = - p_heater_H / des.heater.η_E_H
        elseif des.tes.soc[h,y,s] < controller.options.β_max_tes && h !=1 && !isapprox(des.heater.power_H[h-1,y,s], 0., atol = ϵ)
            p_heater_H = p_net_H
            u_heater_E = - p_heater_H / des.heater.η_E_H
        else
            p_heater_H, u_heater_E = 0., 0.
        end
    else
        p_heater_H, u_heater_E = 0., 0.
    end

    # Store values
    controller.u.liion[h,y,s] = p_net_E - u_heater_E
    controller.u.tes[h,y,s] = p_net_H - p_heater_H
    controller.u.heater[h,y,s] = u_heater_E
    controller.u.elyz[h,y,s] = u_elyz_E
    controller.u.fc[h,y,s] = u_fc_E
    controller.u.h2tank[h,y,s] = u_fc_E / des.fc.η_H2_E + u_elyz_E * des.elyz.η_E_H2
end
function π_2(h::Int64, y::Int64, s::Int64, des::DistributedEnergySystem, controller::RBC)
    controller.u.liion[h,y,s] = des.ld_E.power[h,y,s] - des.pv.power_E[h,y,s]
end
function π_3(h::Int64, y::Int64, s::Int64, des::DistributedEnergySystem, controller::RBC)
    # Control parameters
    β_min_tank, β_max_tank = 0.2, 0.8


    # Net power elec
    p_net_E = des.ld_E.power[h,y,s] - des.pv.power_E[h,y,s]

    if p_net_E < 0
        # Elyz
        if des.h2tank.soc[h,y,s] < β_max_tank
            u_elyz_E, elyz_H, elyz_H2, _ = compute_operation_dynamics(des.elyz, (powerMax = des.elyz.powerMax[y,s], soh = des.elyz.soh[h,y,s]), p_net_E, des.parameters.Δh)
        else
            u_elyz_E, elyz_H, elyz_H2 = 0., 0., 0.
        end
        # Liion
        u_liion = compute_operation_dynamics(des.liion, (Erated = des.liion.Erated[y,s], soc = des.liion.soc[h,y,s], soh = des.liion.soh[h,y,s]), p_net_E - u_elyz_E, des.parameters.Δh)[3]
        # FC
        u_fc_E, fc_H, fc_H2 = 0., 0., 0.
    else
        # Liion
        u_liion = compute_operation_dynamics(des.liion, (Erated = des.liion.Erated[y,s], soc = des.liion.soc[h,y,s], soh = des.liion.soh[h,y,s]), p_net_E, des.parameters.Δh)[3]
        # Elyz
        u_elyz_E, elyz_H, elyz_H2 = 0., 0., 0.
        # FC
        if des.h2tank.soc[h,y,s] > β_min_tank
            u_fc_E, fc_H, fc_H2, _ = compute_operation_dynamics(des.fc, (powerMax = des.fc.powerMax[y,s], soh = des.fc.soh[h,y,s]), p_net_E - u_liion, des.parameters.Δh)
        else
            u_fc_E, fc_H, fc_H2 = 0., 0., 0.
        end
    end

    # Net power heating post H2
    p_net_H = des.ld_H.power[h,y,s] - fc_H - elyz_H

    if p_net_H < 0
        # TODO Rajouter chargement du stockage si surplus elec
        # TES
        u_tes = compute_operation_dynamics(des.tes, (Erated = des.tes.Erated[y], soc = des.tes.soc[h,y,s]), p_net_H, des.parameters.Δh)[2]
        u_heater_E, heater_H = 0., 0.
    else
        # TODO Si surplus elec, alors directement les utiliser pour couvrir le besoin en chaud !
        # TES
        u_tes = compute_operation_dynamics(des.tes, (Erated = des.tes.Erated[y], soc = des.tes.soc[h,y,s]), p_net_H, des.parameters.Δh)[2]
        # Heater
        u_heater_E, heater_H = compute_operation_dynamics(des.heater, (powerMax = des.heater.powerMax[y,s],), - (p_net_H - u_tes) / des.heater.η_E_H, des.parameters.Δh)
    end

    # Store values
    controller.u.liion[h,y,s] = p_net_E - u_fc_E - u_elyz_E - u_heater_E
    controller.u.tes[h,y,s] = u_tes
    controller.u.heater[h,y,s] = u_heater_E
    controller.u.elyz[h,y,s] = u_elyz_E
    controller.u.fc[h,y,s] = u_fc_E
    controller.u.h2tank[h,y,s] = compute_operation_dynamics(des.h2tank, (Erated = des.h2tank.Erated[y,s], soc = des.h2tank.soc[h,y,s]), - fc_H2 - elyz_H2, des.parameters.Δh)[2]
end

### Offline
function initialize_controller!(des::DistributedEnergySystem, controller::RBC, ω::AbstractScenarios)
    # Save history for online optimization
    # controller.history = ω

    # Preallocation
    preallocate!(controller, des.parameters.nh, des.parameters.ny, des.parameters.ns)

    return controller
end

### Online
function compute_operation_decisions!(h::Int64, y::Int64, s::Int64, des::DistributedEnergySystem, controller::RBC)
    # Chose policy TODO : better way !
    if isa(des.ld_H, Load)
        return π_3(h, y, s, des, controller)
    else
        return π_2(h, y, s, des, controller)
    end
end
