#=
    This file includes all the funtions needed to compute the operation
    and investment dynamics
 =#

function compute_operation_dynamics!(h::Int64, y::Int64, s::Int64, des::DistributedEnergySystem, controller::AbstractController)

    # Converters
    if isa(des.heater, Heater)
        des.heater.power_E[h,y,s], des.heater.power_H[h,y,s] =
        compute_operation_dynamics(des.heater, (powerMax = des.heater.powerMax[y,s],), controller.u.heater[h,y,s], des.parameters.Δh)
    end

    if isa(des.elyz, Electrolyzer)
        des.elyz.power_E[h,y,s], des.elyz.power_H[h,y,s], des.elyz.power_H2[h,y,s], des.elyz.soh[h+1,y,s] =
        compute_operation_dynamics(des.elyz, (powerMax = des.elyz.powerMax[y,s], soh = des.elyz.soh[h,y,s]), controller.u.elyz[h,y,s], des.parameters.Δh)
    end

    if isa(des.fc, FuelCell)
        des.fc.power_E[h,y,s], des.fc.power_H[h,y,s], des.fc.power_H2[h,y,s], des.fc.soh[h+1,y,s] =
        compute_operation_dynamics(des.fc, (powerMax = des.fc.powerMax[y,s], soh = des.fc.soh[h,y,s]), controller.u.fc[h,y,s], des.parameters.Δh)
    end

    # Storage
    if isa(des.liion, Liion)
        des.liion.soc[h+1,y,s], des.liion.soh[h+1,y,s], des.liion.power_E[h,y,s] =
        compute_operation_dynamics(des.liion, (Erated = des.liion.Erated[y,s], soc = des.liion.soc[h,y,s], soh = des.liion.soh[h,y,s]), controller.u.liion[h,y,s], des.parameters.Δh)
    end

    if isa(des.tes, ThermalSto)
        des.tes.soc[h+1,y,s], des.tes.power_H[h,y,s] =
        compute_operation_dynamics(des.tes, (Erated = des.tes.Erated[y], soc = des.tes.soc[h,y,s]), controller.u.tes[h,y,s], des.parameters.Δh)
    end

    if isa(des.h2tank, H2Tank)
        des.h2tank.soc[h+1,y,s], des.h2tank.power_H2[h,y,s] =
        compute_operation_dynamics(des.h2tank, (Erated = des.h2tank.Erated[y,s], soc = des.h2tank.soc[h,y,s]), controller.u.h2tank[h,y,s], des.parameters.Δh)
    end
end
function compute_investment_dynamics!(y::Int64, s::Int64, des::DistributedEnergySystem, designer::AbstractDesigner)

     # Converters
     if isa(des.pv, Source)
         des.pv.powerMax[y+1,s] =
         compute_investment_dynamics(des.pv, (powerMax = des.pv.powerMax[y,s],), designer.u.pv[y,s])
     end

     if isa(des.elyz, Electrolyzer)
         des.elyz.powerMax[y+1,s], des.elyz.soh[1,y+1,s] =
         compute_investment_dynamics(des.elyz, (powerMax = des.elyz.powerMax[y,s], soh = des.elyz.soh[end,y,s]), designer.u.elyz[y,s])
     end

     if isa(des.fc, FuelCell)
         des.fc.powerMax[y+1,s], des.fc.soh[1,y+1,s] =
         compute_investment_dynamics(des.fc, (powerMax = des.fc.powerMax[y,s], soh = des.fc.soh[end,y,s]), designer.u.fc[y,s])
     end

     # Storage
     if isa(des.liion, Liion)
         des.liion.Erated[y+1,s], des.liion.soc[1,y+1,s], des.liion.soh[1,y+1,s] =
         compute_investment_dynamics(des.liion, (Erated = des.liion.Erated[y,s], soc = des.liion.soc[end,y,s], soh = des.liion.soh[end,y,s]), designer.u.liion[y,s])
     end

     if isa(des.tes, ThermalSto)
         des.tes.Erated[y+1,s], des.tes.soc[1,y+1,s] =
         compute_investment_dynamics(des.tes, (Erated = des.tes.Erated[y,s], soc = des.tes.soc[end,y,s]), designer.u.tes[y,s])
     end

     if isa(des.h2tank, H2Tank)
         des.h2tank.Erated[y+1,s], des.h2tank.soc[1,y+1,s] =
         compute_investment_dynamics(des.h2tank, (Erated = des.h2tank.Erated[y,s], soc = des.h2tank.soc[end,y,s]), designer.u.h2tank[y,s])
     end
end
