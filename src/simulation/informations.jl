#=
    This file includes all the funtions needed to update the operation
    and investment informations
 =#


function update_operation_informations!(h::Int64, y::Int64, s::Int64, des::DES, ω::AbstractScenarios)

    isa(des.ld_E, Load) ? des.ld_E.power[h,y,s] = ω.values.ld_E[h,y,s] : nothing

    isa(des.ld_H, Load) ? des.ld_H.power[h,y,s] = ω.values.ld_H[h,y,s] : nothing

    isa(des.pv, Source) ? des.pv.power_E[h,y,s] = des.pv.powerMax[y,s] * ω.values.pv_E[h,y,s] : nothing

    if isa(des.grid, Grid)
        des.grid.C_grid_in[h,y,s] = ω.values.C_grid_in[h,y,s]
        des.grid.C_grid_out[h,y,s] = ω.values.C_grid_out[h,y,s]
    end
end
function update_investment_informations!(y::Int64, s::Int64, des::DES, ω::AbstractScenarios)

    isa(des.pv, Source) ? des.pv.C_pv[y,s] = ω.values.C_pv[y,s] : nothing

    isa(des.liion, Liion) ? des.liion.C_liion[y,s] = ω.values.C_liion[y,s] : nothing

    isa(des.tes, ThermalSto) ? des.tes.C_tes[y,s] = ω.values.C_tes[y,s] : nothing

    isa(des.h2tank, H2Tank) ? des.h2tank.C_tank[y,s] = ω.values.C_tank[y,s] : nothing

    isa(des.elyz, Electrolyzer) ? des.elyz.C_elyz[y,s] = ω.values.C_elyz[y,s] : nothing

    isa(des.fc, FuelCell) ? des.fc.C_fc[y,s] = ω.values.C_fc[y,s] : nothing

    isa(des.heater, Heater) ? des.heater.C_heater[y,s] = ω.values.C_heater[y,s] : nothing
end
