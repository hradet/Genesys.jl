#=
    Rule based controller
=#

mutable struct RuleBasedController <: AbstractController
    u::NamedTuple
    π::Function
    parameters::Dict{String, Any}
    RuleBasedController() = new()
end

#### Policies definition ####
# Simple
function rb_operation_policy(h::Int64, y::Int64, s::Int64, ld::Load, pv::Source, liion::Liion, controller::RuleBasedController)
    # Liion
    u_liion = ld.power[h,y,s] - pv.power_E[h,y,s]
    return u_liion
end
# Multi-energy
function rb_operation_policy(h::Int64, y::Int64, s::Int64, ld::Load, pv::Source, liion::Liion,
    h2tank::H2Tank, elyz::Electrolyzer, fc::FuelCell, tes::ThermalSto, heater::Heater, controller::RuleBasedController)
    # Control parameters
    ϵ = 1e-2

    # Net power elec
    p_net_E = ld.power_E[h,y,s] - pv.power_E[h,y,s]

    # H2 - Ullberg, 2004
    if p_net_E <= 0
        # FC is off
        u_fc, p_fc_H = 0., 0.
        # Strategy for Elyz
        if liion.soc[h,y,s] > controller.parameters["β_max_elyz"]
            u_elyz = min(max(p_net_E, -elyz.powerMax[y,s]), -elyz.α_p * elyz.powerMax[y,s])
            p_elyz_H = - u_elyz * elyz.η_E_H
        elseif liion.soc[h,y,s] > controller.parameters["β_min_elyz"] && h !=1 && elyz.power_E[h-1,y,s] <= -ϵ
            u_elyz = min(max(p_net_E, -elyz.powerMax[y,s]), -elyz.α_p * elyz.powerMax[y,s])
            p_elyz_H = - u_elyz * elyz.η_E_H
        else
            u_elyz, p_elyz_H = 0., 0.
        end
    else
        # Elyz is off
        u_elyz, p_elyz_H = 0., 0.
        # Strategy for FC
        if liion.soc[h,y,s] < controller.parameters["β_min_fc"]
            u_fc = min(max(p_net_E, fc.α_p * fc.powerMax[y,s]), fc.powerMax[y,s])
            p_fc_H = u_fc * fc.η_H2_H / fc.η_H2_E
        elseif liion.soc[h,y,s] < controller.parameters["β_max_fc"] && h !=1 && fc.power_E[h-1,y,s] >= ϵ
            u_fc = min(max(p_net_E, fc.α_p * fc.powerMax[y,s]), fc.powerMax[y,s])
            p_fc_H = u_fc * fc.η_H2_H / fc.η_H2_E
        else
            u_fc, p_fc_H = 0., 0.
        end
    end

    # Net power elec post H2
    p_net_E += - u_elyz - u_fc

    # Net power heating post H2
    p_net_H = ld.power_H[h,y,s] - p_fc_H - p_elyz_H

    # Heater
    if p_net_H >= 0.
        if tes.soc[h,y,s] < controller.parameters["β_min_tes"]
            p_heater_H = p_net_H
            u_heater = - p_heater_H / heater.η_E_H
        elseif tes.soc[h,y,s] < controller.parameters["β_max_tes"] && h !=1 && heater.power_H[h-1,y,s] > ϵ
            p_heater_H = p_net_H
            u_heater = - p_heater_H / heater.η_E_H
        else
            p_heater_H, u_heater = 0., 0.
        end
    else
        p_heater_H, u_heater = 0., 0.
    end

    # Liion
    u_liion = p_net_E - u_heater

    # TES
    u_tes = p_net_H - p_heater_H

    return u_liion, u_elyz, u_fc, u_tes, u_heater
end

#### Offline functions ####
# Simple
function initialize_controller(ld::Load, pv::Source, liion::Liion,
    controller::RuleBasedController, grid::Grid, ω_optim::Scenarios, parameters::NamedTuple)
     # Parameters
     nh = size(ld.power_E,1) # number of simulation hours in one year
     ny = size(ld.power_E,2) # number of simulation years
     ns = size(ld.power_E,3) # number of scenarios

     # Initialize controller policy
     controller.π = rb_operation_policy

     # Initialize decisions variables
     controller.u = (
     u_liion = convert(SharedArray,zeros(nh,ny,ns)),
     )
end
# Multi-energy
function initialize_controller(ld::Load, pv::Source, liion::Liion, h2tank::H2Tank,
      elyz::Electrolyzer, fc::FuelCell, tes::ThermalSto, heater::Heater,
      controller::RuleBasedController, grid::Grid, ω_optim::Scenarios, parameters::NamedTuple)
      # Parameters
      nh = size(ld.power_E,1) # number of simulation hours in one year
      ny = size(ld.power_E,2) # number of simulation years
      ns = size(ld.power_E,3) # number of scenarios

      # Initialize controller policy
      controller.π = rb_operation_policy

      # Initialize decisions variables
      # Operation decisions
      controller.u = (
      u_liion = zeros(nh,ny,ns),
      u_elyz = zeros(nh,ny,ns),
      u_fc = zeros(nh,ny,ns),
      u_tes = zeros(nh,ny,ns),
      u_heater = zeros(nh,ny,ns),
      )
end

#### Online functions ####
# Simple
function compute_operation_decisions(h::Int64, y::Int64, s::Int64, ld::Load, pv::Source,
     liion::Liion, grid::Grid, controller::RuleBasedController, ω_optim::Scenarios, parameters::NamedTuple)
    controller.u.u_liion[h,y,s] = controller.π(h, y, s, ld, pv, liion, controller)
end
# Multi-energy
function compute_operation_decisions(h::Int64, y::Int64, s::Int64, ld::Load, pv::Source,
     liion::Liion, h2tank::H2Tank, elyz::Electrolyzer, fc::FuelCell, tes::ThermalSto,
     heater::Heater, controller::RuleBasedController, ω_optim::Scenarios, parameters::NamedTuple)
    controller.u.u_liion[h,y,s], controller.u.u_elyz[h,y,s], controller.u.u_fc[h,y,s], controller.u.u_tes[h,y,s], controller.u.u_heater[h,y,s] = controller.π(h, y, s, ld, pv, liion, h2tank, elyz, fc, tes, heater, controller)
end
