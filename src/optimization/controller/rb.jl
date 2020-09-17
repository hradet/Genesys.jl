#=
    Rule based controller
=#
mutable struct RBCOptions
    β_min_tes
    β_max_tes
    β_min_fc
    β_max_fc
    β_min_elyz
    β_max_elyz

    RBCOptions(; β_min_tes = 0.2,
                       β_max_tes = 0.9,
                       β_min_fc = 0.25,
                       β_max_fc = 0.3,
                       β_min_elyz = 0.4,
                       β_max_elyz = 0.45) =
                       new(β_min_tes, β_max_tes, β_min_fc, β_max_fc, β_min_elyz, β_max_elyz)
end

mutable struct RBC <: AbstractController
    options::RBCOptions
    u::NamedTuple
    π::Function
    RBC(; options = RBCOptions()) = new(options)
end

### Policies
function policy_1(h::Int64, y::Int64, s::Int64, des::DistributedEnergySystem, controller::RBC)
    # Control parameters
    ϵ = 1e-2

    # Net power elec
    p_net_E = des.ld_E.power[h,y,s] - des.pv.power_E[h,y,s]

    # H2 - Ullberg, 2004
    if p_net_E <= 0
        # FC is off
        u_fc_E, p_fc_H = 0., 0.
        # Strategy for Elyz
        if des.liion.soc[h,y,s] > controller.options.β_max_elyz
            u_elyz_E = min(max(p_net_E, -des.elyz.powerMax[y,s]), -des.elyz.α_p * des.elyz.powerMax[y,s])
            p_elyz_H = - u_elyz_E * des.elyz.η_E_H
        elseif des.liion.soc[h,y,s] > controller.options.β_min_elyz && h !=1 && des.elyz.power_E[h-1,y,s] <= -ϵ
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
        elseif des.liion.soc[h,y,s] < controller.options.β_max_fc && h !=1 && des.fc.power_E[h-1,y,s] >= ϵ
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
        elseif des.tes.soc[h,y,s] < controller.options.β_max_tes && h !=1 && des.heater.power_H[h-1,y,s] > ϵ
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
function policy_2(h::Int64, y::Int64, s::Int64, des::DistributedEnergySystem, controller::RBC)
    controller.u.liion[h,y,s] = des.ld_E.power[h,y,s] - des.pv.power_E[h,y,s]
end

### Offline
function initialize_controller!(des::DistributedEnergySystem, controller::RBC, ω::AbstractScenarios)
    # Preallocation
    preallocate!(controller, des.parameters.nh, des.parameters.ny, des.parameters.ns)

    # Chose policy TODO : better way !
    if isa(des.ld_H, Load)
        controller.π = policy_1
    else
        controller.π = policy_2
    end

    return controller
end

### Online
function compute_operation_decisions!(h::Int64, y::Int64, s::Int64, des::DistributedEnergySystem, controller::RBC)
    controller.π(h, y, s, des, controller)
end
