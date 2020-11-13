#=
    This file includes all the funtions needed to update the operation
    and investment informations
 =#


function update_operation_informations!(h::Int64, y::Int64, s::Int64, des::DistributedEnergySystem, ω::AbstractScenarios)

    if isa(des.ld_E, Load)
        des.ld_E.power[h,y,s] = ω.ld_E.power[h,y,s]
        des.ld_E.timestamp[h,y,s] = ω.ld_E.t[h,y,s]
    end

    if isa(des.ld_H, Load)
         des.ld_H.power[h,y,s] = ω.ld_H.power[h,y,s]
         des.ld_H.timestamp[h,y,s] = ω.ld_E.t[h,y,s]
     end

    if isa(des.pv, Source)
        des.pv.power_E[h,y,s] = des.pv.powerMax[y,s] * ω.pv.power[h,y,s]
        des.pv.timestamp[h,y,s] = ω.pv.t[h,y,s]
    end

    if isa(des.grid, Grid)
        des.grid.cost_in[h,y,s] = ω.grid.cost_in[h,y,s]
        des.grid.cost_out[h,y,s] = ω.grid.cost_out[h,y,s]
    end
end
function update_investment_informations!(y::Int64, s::Int64, des::DistributedEnergySystem, ω::AbstractScenarios)

    isa(des.pv, Source) ? des.pv.cost[y,s] = ω.pv.cost[y,s] : nothing

    isa(des.liion, Liion) ? des.liion.cost[y,s] = ω.liion.cost[y,s] : nothing

    isa(des.tes, ThermalSto) ? des.tes.cost[y,s] = ω.tes.cost[y,s] : nothing

    isa(des.h2tank, H2Tank) ? des.h2tank.cost[y,s] = ω.h2tank.cost[y,s] : nothing

    isa(des.elyz, Electrolyzer) ? des.elyz.cost[y,s] = ω.elyz.cost[y,s] : nothing

    isa(des.fc, FuelCell) ? des.fc.cost[y,s] = ω.fc.cost[y,s] : nothing

    isa(des.heater, Heater) ? des.heater.cost[y,s] = ω.heater.cost[y,s] : nothing
end
