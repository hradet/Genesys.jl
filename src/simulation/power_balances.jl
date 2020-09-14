#=
    This file includes all the functions needed to check the power balance
    constraints
=#

function compute_power_balances!(h::Int64, y::Int64, s::Int64, des::DES)

    # Hydrogen
    check_hydrogen!(h, y, s, des)

    # Heat
    check_heat!(h, y, s, des)

    # Electricity
    check_electricity!(h, y, s, des)
end
function check_hydrogen!(h::Int64, y::Int64, s::Int64, des::DES)

    isa(des.h2tank, H2Tank) ? h2tank = des.h2tank.power_H2[h,y,s] : h2tank = 0.
    isa(des.elyz, Electrolyzer) ? elyz = des.elyz.power_H2[h,y,s] : elyz = 0.
    isa(des.fc, FuelCell) ? fc = des.fc.power_H2[h,y,s] : fc = 0.

    if h2tank + elyz +  fc < 0.

        if isa(des.h2tank, H2Tank)
            des.h2tank.power_H2[h,y,s] = 0.
            des.h2tank.soc[h+1,y,s] = max(0., des.h2tank.soc[h,y,s] * (1. - des.h2tank.η_self))
        end

        if isa(des.elyz, Electrolyzer)
            des.elyz.power_E[h,y,s], des.elyz.power_H[h,y,s], des.elyz.power_H2[h,y,s]  = 0., 0., 0.
            des.elyz.soh[h+1,y,s] = des.elyz.soh[h,y,s]
        end

        if isa(des.fc, FuelCell)
            des.fc.power_E[h,y,s], des.fc.power_H[h,y,s], des.fc.power_H2[h,y,s] = 0., 0., 0.
            des.fc.soh[h+1,y,s] = des.fc.soh[h,y,s]
        end

    end
end
function check_heat!(h::Int64, y::Int64, s::Int64, des::DES)

    isa(des.ld_H, Load) ? ld_H = des.ld_H.power_H[h,y,s] : ld_H = 0.
    isa(des.heater, Heater) ? heater = des.heater.power_H[h,y,s] : heater = 0.
    isa(des.tes, ThermalSto) ? tes = des.tes.power_H[h,y,s] : tes = 0.
    isa(des.elyz, Electrolyzer) ? elyz = des.elyz.power_H[h,y,s] : elyz = 0.
    isa(des.fc, FuelCell) ? fc = des.fc.power_H[h,y,s] : fc = 0.

    if heater + tes + elyz +  fc < ld_H
        if heater + tes + elyz +  fc < 0.
            if isa(des.tes, ThermalSto)
                des.tes.power_H[h,y,s] = 0.
                des.tes.soc[h+1,y,s] = max(0., des.tes.soc[h,y,s] * (1. - des.tes.η_self))
            end

            if isa(des.heater, Heater)
                des.heater.power_E[h,y,s], des.heater.power_H[h,y,s] = 0., 0.
            end
        end
    else
        # print("Warning! Thermal load not fully supplied")
    end
end
function check_electricity!(h::Int64, y::Int64, s::Int64, des::DES)

    isa(des.ld_E, Load) ? ld_E = des.ld_E.power[h,y,s] : ld_E = 0.
    isa(des.pv, Source) ? pv = des.pv.power_E[h,y,s] : pv = 0.
    isa(des.heater, Heater) ? heater = des.heater.power_E[h,y,s] : heater = 0.
    isa(des.liion, Liion) ? liion = des.liion.power_E[h,y,s] : liion = 0.
    isa(des.elyz, Electrolyzer) ? elyz = des.elyz.power_E[h,y,s] : elyz = 0.
    isa(des.fc, FuelCell) ? fc = des.fc.power_E[h,y,s] : fc = 0.
    # Power grid
    des.grid.power_E[h,y,s] = ld_E - pv - liion - elyz - fc - heater
end
