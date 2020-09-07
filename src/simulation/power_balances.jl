#=
    This file includes all the functions needed to check the power balance
    constraints
=#
# Simple
function power_balance(h::Int64, y::Int64, s::Int64, ld::Load, pv::Source, liion::Liion, grid::Grid)
    # Electric power balance
    grid.power_E[h,y,s] = max(0. , ld.power_E[h,y,s] - pv.power_E[h,y,s] - liion.power_E[h,y,s])
end
# Multi-energy
# TODO : diviser les powers balance pour chaque noeud d'énergie + in/out réseau
function power_balance(h::Int64, y::Int64, s::Int64, ld::Load, pv::Source, liion::Liion,
    h2tank::H2Tank, elyz::Electrolyzer, fc::FuelCell, tes::ThermalSto, heater::Heater, grid::Grid)

    # Hydrogen power balance
    if h2tank.power_H2[h,y,s] + elyz.power_H2[h,y,s] +  fc.power_H2[h,y,s] < 0.

        # Power correction
        h2tank.power_H2[h,y,s] = 0.
        elyz.power_E[h,y,s], elyz.power_H[h,y,s], elyz.power_H2[h,y,s]  = 0., 0., 0.
        fc.power_E[h,y,s], fc.power_H[h,y,s], fc.power_H2[h,y,s] = 0., 0., 0.

        # SoC and SoH correction
        h2tank.soc[h+1,y,s] = max(0., h2tank.soc[h,y,s] * (1. - h2tank.η_self))
        elyz.soh[h+1,y,s] = elyz.soh[h,y,s]
        fc.soh[h+1,y,s] = fc.soh[h,y,s]

    end

    # Thermal power balance
    if heater.power_H[h,y,s] + tes.power_H[h,y,s] + elyz.power_H[h,y,s] +  fc.power_H[h,y,s] < ld.power_H[h,y,s]

        if heater.power_H[h,y,s] + tes.power_H[h,y,s] + elyz.power_H[h,y,s] +  fc.power_H[h,y,s] < 0.

            # Power correction
            heater.power_H[h,y,s], heater.power_E[h,y,s]  = 0., 0.
            tes.power_H[h,y,s] = 0.

            # SoC and SoH correction
            tes.soc[h+1,y,s] = max(0., tes.soc[h,y,s] * (1. - tes.η_self))

        else

            # print("Warning! Thermal load not fully supplied")

        end
    end

    # Electric power balance
    grid.power_E[h,y,s] = max(0. , ld.power_E[h,y,s] - pv.power_E[h,y,s] - liion.power_E[h,y,s] - elyz.power_E[h,y,s] - fc.power_E[h,y,s] - heater.power_E[h,y,s])
end
