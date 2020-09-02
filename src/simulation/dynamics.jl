#=
    This file includes all the funtions needed to compute the operation
    and investment dynamics
 =#

# Simple
function compute_operation_dynamics(h::Int64, y::Int64, s::Int64, liion::Liion, controller::AbstractController, parameters::NamedTuple)
    liion.soc[h+1,y,s], liion.soh[h+1,y,s], liion.power_E[h,y,s] = compute_operation_dynamics(liion, (Erated = liion.Erated[y,s], soc = liion.soc[h,y,s], soh = liion.soh[h,y,s]), controller.u.u_liion[h,y,s], parameters.Δh)
end
function compute_investment_dynamics(y::Int64, s::Int64, pv::Source, liion::Liion, designer::AbstractDesigner)
     # Converters
     pv.powerMax[y+1,s] = compute_investment_dynamics(pv, (powerMax = pv.powerMax[y,s],), designer.u.u_pv[y,s])
     # Storage
     liion.Erated[y+1,s], liion.soc[1,y+1,s], liion.soh[1,y+1,s] = compute_investment_dynamics(liion, (Erated = liion.Erated[y,s], soc = liion.soc[end,y,s], soh = liion.soh[end,y,s]), designer.u.u_liion[y,s])
end
# Multi-energy
function compute_operation_dynamics(h::Int64, y::Int64, s::Int64, liion::Liion, h2tank::H2Tank,
    elyz::Electrolyzer, fc::FuelCell, tes::ThermalSto, heater::Heater, controller::AbstractController, parameters::NamedTuple)

    # Converters
    heater.power_E[h,y,s], heater.power_H[h,y,s] = compute_operation_dynamics(heater, (powerMax = heater.powerMax[y,s],), controller.u.u_heater[h,y,s], parameters.Δh)
    elyz.power_E[h,y,s], elyz.power_H[h,y,s], elyz.power_H2[h,y,s], elyz.soh[h+1,y,s] = compute_operation_dynamics(elyz, (powerMax = elyz.powerMax[y,s], soh = elyz.soh[h,y,s]), controller.u.u_elyz[h,y,s], parameters.Δh)
    fc.power_E[h,y,s], fc.power_H[h,y,s], fc.power_H2[h,y,s], fc.soh[h+1,y,s] = compute_operation_dynamics(fc, (powerMax = fc.powerMax[y,s], soh = fc.soh[h,y,s]), controller.u.u_fc[h,y,s], parameters.Δh)

    # Storage
    liion.soc[h+1,y,s], liion.soh[h+1,y,s], liion.power_E[h,y,s] = compute_operation_dynamics(liion, (Erated = liion.Erated[y,s], soc = liion.soc[h,y,s], soh = liion.soh[h,y,s]), controller.u.u_liion[h,y,s], parameters.Δh)
    tes.soc[h+1,y,s], tes.power_H[h,y,s] = compute_operation_dynamics(tes, (Erated = tes.Erated[y], soc = tes.soc[h,y,s]), controller.u.u_tes[h,y,s], parameters.Δh)
    h2tank.soc[h+1,y,s], h2tank.power_H2[h,y,s] = compute_operation_dynamics(h2tank, (Erated = h2tank.Erated[y,s], soc = h2tank.soc[h,y,s]), - (fc.power_H2[h,y,s] + elyz.power_H2[h,y,s]), parameters.Δh)
end
function compute_investment_dynamics(y::Int64, s::Int64, pv::Source, liion::Liion, h2tank::H2Tank,
     elyz::Electrolyzer, fc::FuelCell, tes::ThermalSto, designer::AbstractDesigner)

     # Converters
     pv.powerMax[y+1,s] = compute_investment_dynamics(pv, (powerMax = pv.powerMax[y,s],), designer.u.u_pv[y,s])
     elyz.powerMax[y+1,s], elyz.soh[1,y+1,s] = compute_investment_dynamics(elyz, (powerMax = elyz.powerMax[y,s], soh = elyz.soh[end,y,s]), designer.u.u_elyz[y,s])
     fc.powerMax[y+1,s], fc.soh[1,y+1,s] = compute_investment_dynamics(fc, (powerMax = fc.powerMax[y,s], soh = fc.soh[end,y,s]), designer.u.u_fc[y,s])
     # Storage
     liion.Erated[y+1,s], liion.soc[1,y+1,s], liion.soh[1,y+1,s] = compute_investment_dynamics(liion, (Erated = liion.Erated[y,s], soc = liion.soc[end,y,s], soh = liion.soh[end,y,s]), designer.u.u_liion[y,s])
     tes.Erated[y+1,s], tes.soc[1,y+1,s] = compute_investment_dynamics(tes, (Erated = tes.Erated[y,s], soc = tes.soc[end,y,s]), designer.u.u_tes[y,s])
     h2tank.Erated[y+1,s], h2tank.soc[1,y+1,s] = compute_investment_dynamics(h2tank, (Erated = h2tank.Erated[y,s], soc = h2tank.soc[end,y,s]), designer.u.u_tank[y,s])
end
