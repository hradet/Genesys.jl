#=
    This file includes all the functions needed to check the power balance
    constraints
=#

function power_balance_checking(h, y, s, ld, pv, liion, h2tank, elyz, fc, tes, heater, grid)

    # Hydrogen power balance
    if h2tank.power_H2[h,y,s] + elyz.power_H2[h,y,s] +  fc.power_H2[h,y,s] < 0.

        # Power correction
        h2tank.power_H2[h,y,s] = 0.
        elyz.power_E[h,y,s] = 0.
        elyz.power_H[h,y,s] = 0.
        elyz.power_H2[h,y,s] = 0.
        fc.power_E[h,y,s] = 0.
        fc.power_H[h,y,s] = 0.
        fc.power_H2[h,y,s] = 0.

        # SoC and SoH correction
        h2tank.soc[h+1,y,s] = max(0., h2tank.soc[h,y,s] * (1 - h2tank.η_self))
        elyz.soh[h+1,y,s] = elyz.soh[h,y,s]
        fc.soh[h+1,y,s] = fc.soh[h,y,s]

    end

    # Thermal power balance
    if heater.power_H[h,y,s] + tes.power_H[h,y,s] + elyz.power_H[h,y,s] +  fc.power_H[h,y,s] < ld.power_H[h,y,s]

        if heater.power_H[h,y,s] + tes.power_H[h,y,s] + elyz.power_H[h,y,s] +  fc.power_H[h,y,s] < 0.

            # Power correction
            heater.power_H[h,y,s] = 0.
            heater.power_E[h,y,s]  = 0.
            tes.power_H[h,y,s] = 0.

            # SoC and SoH correction
            tes.soc[h+1,y,s] = max(0., tes.soc[h,y,s] * (1 - tes.η_self))

        else

            # print("Warning! Thermal load not fully supplied")

        end
    end

    # Electric power balance
    grid.power_E[h,y,s] = max(0. , ld.power_E[h,y,s] - pv.power_E[h,y,s] - liion.power_E[h,y,s] - elyz.power_E[h,y,s] - fc.power_E[h,y,s] - heater.power_E[h,y,s])
end

# TODO : diviser les powers balance pour chaque noeud d'énergie
