#=
    This file includes all the functions needed to check the power balance
    constraints
=#

#TODO rajouter fonction add_to_powerbalance!(power_balance, power)

function compute_power_balances!(h::Int64, y::Int64, s::Int64, des::DistributedEnergySystem)

    # Hydrogen
    check_hydrogen!(h, y, s, des)

    # Heat
    check_heat!(h, y, s, des)

    # Electricity
    check_electricity!(h, y, s, des)
end
function check_hydrogen!(h::Int64, y::Int64, s::Int64, des::DistributedEnergySystem)

    ϵ = 0.01 # 0.01 kW tolerance

    isa(des.h2tank, H2Tank) ? h2tank = des.h2tank.power_H2[h,y,s] : h2tank = 0.
    isa(des.elyz, Electrolyzer) ? elyz = des.elyz.power_H2[h,y,s] : elyz = 0.
    isa(des.fc, FuelCell) ? fc = des.fc.power_H2[h,y,s] : fc = 0.

    if !isapprox(h2tank + elyz +  fc, 0., atol = ϵ)

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
function check_heat!(h::Int64, y::Int64, s::Int64, des::DistributedEnergySystem)

    ϵ = 0.01 # 0.01 kW tolerance

    isa(des.ld_H, Load) ? ld_H = des.ld_H.power[h,y,s] : ld_H = 0.
    isa(des.heater, Heater) ? heater = des.heater.power_H[h,y,s] : heater = 0.
    isa(des.tes, ThermalSto) ? tes = des.tes.power_H[h,y,s] : tes = 0.
    isa(des.elyz, Electrolyzer) ? elyz = des.elyz.power_H[h,y,s] : elyz = 0.
    isa(des.fc, FuelCell) ? fc = des.fc.power_H[h,y,s] : fc = 0.

    if !isapprox(heater + tes + elyz +  fc, ld_H, atol = ϵ)
        if !isapprox(heater + tes + elyz +  fc, 0., atol = ϵ)
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
function check_electricity!(h::Int64, y::Int64, s::Int64, des::DistributedEnergySystem)

    isa(des.ld_E, Load) ? ld_E = des.ld_E.power[h,y,s] : ld_E = 0.
    isa(des.pv, Source) ? pv = des.pv.power_E[h,y,s] : pv = 0.
    isa(des.heater, Heater) ? heater = des.heater.power_E[h,y,s] : heater = 0.
    isa(des.liion, Liion) ? liion = des.liion.power_E[h,y,s] : liion = 0.
    isa(des.elyz, Electrolyzer) ? elyz = des.elyz.power_E[h,y,s] : elyz = 0.
    isa(des.fc, FuelCell) ? fc = des.fc.power_E[h,y,s] : fc = 0.
    # Power grid
    isa(des.grid, Grid) ? des.grid.power_E[h,y,s] = max(min(ld_E - pv - liion - elyz - fc - heater, des.grid.powerMax), -des.grid.powerMax)  : nothing
end
