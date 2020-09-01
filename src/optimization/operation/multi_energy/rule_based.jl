#=
    Rule based controller
=#

#                                   Policy definition
#_______________________________________________________________________________
function rb_operation_policy(h::Int64, y::Int64, s::Int64, ld::Load, pv::Source, liion::Liion,
    h2tank::H2Tank, elyz::Electrolyzer, fc::FuelCell, tes::ThermalSto, heater::Heater)
    # Control parameters
    ϵ = 1e-2
    threshold_min_tes, threshold_max_tes = 0.1, 0.9
    threshold_min_h2, threshold_max_h2 = 0.2, 0.8
    threshold_min_liion, threshold_max_liion = 0.3, 0.7

    # Heater
    if tes.soc[h,y,s] < threshold_min_tes
        powerHeater_H = ld.power_H[h,y,s]
        u_heater = - powerHeater_H / heater.η_E_H
    elseif tes.soc[h,y,s] < threshold_max_tes && h !=1 && heater.power_H[h-1,y,s] > ϵ
        powerHeater_H = ld.power_H[h,y,s]
        u_heater = - powerHeater_H / heater.η_E_H
    else
        powerHeater_H = u_heater = 0.
    end

    # Net power elec
    powerNet_E = ld.power_E[h,y,s] - pv.powerMax[y,s] * pv.power_E[h,y,s] - u_heater

    # H2Hub
    if powerNet_E < -elyz.α_p * elyz.powerMax[y,s]
        u_fc = 0
        if h2tank.soc[h,y,s] < threshold_max_h2
            if liion.soc[h,y,s] > threshold_max_liion
                u_elyz = powerNet_E
                powerH2_H = - u_elyz * elyz.η_E_H
            elseif liion.soc[h,y,s] > threshold_min_liion && h !=1 && elyz.power_E[h-1,y,s] <= -ϵ
                u_h2hub = powerNet_E
                powerH2_H = - u_elyz * elyz.η_E_H
            else
                u_elyz = powerH2_H = 0.
            end
        else
            u_elyz = powerH2_H = 0.
        end
    elseif powerNet_E > fc.α_p * fc.powerMax[y,s]
        u_elyz = 0
        if h2tank.soc[h,y,s] > threshold_min_h2
            u_fc = powerNet_E
            powerH2_H = u_fc * fc.η_H2_H / fc.η_H2_E
        else
            u_fc = powerH2_H = 0.
        end
    else
        u_fc = u_elyz = powerH2_H = 0.
    end

    # Net power elec post H2
    powerNet_E = powerNet_E - u_elyz - u_fc

    # Liion
    u_liion = powerNet_E

    # Net power heating post H2
    powerNet_H = ld.power_H[h,y,s] - powerH2_H - powerHeater_H

    # TES
    u_tes = powerNet_H

    return u_liion, u_elyz, u_fc, u_tes, u_heater
end

#                                   Offline functions
#_______________________________________________________________________________
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

#                                   Online functions
#_______________________________________________________________________________
function compute_operation_decisions(h::Int64, y::Int64, s::Int64, ld::Load, pv::Source,
     liion::Liion, h2tank::H2Tank, elyz::Electrolyzer, fc::FuelCell, tes::ThermalSto,
     heater::Heater, controller::RuleBasedController, ω_optim::Scenarios, parameters::NamedTuple)
    controller.u.u_liion[h,y,s], controller.u.u_elyz[h,y,s], controller.u.u_fc[h,y,s], controller.u.u_tes[h,y,s], controller.u.u_heater[h,y,s] = controller.π(h, y, s, ld, pv, liion, h2tank, elyz, fc, tes, heater)
end
