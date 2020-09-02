#=
    This file includes all the funtions needed to update the operation
    and investment informations
 =#

# Simple
function update_operation_informations(h::Int64, y::Int64, s::Int64, ld::Load, pv::Source, liion::Liion, grid::Grid, ω::AbstractScenarios)
    # Energy demand and production
    ld.power_E[h,y,s] = ω.ld_E[h,y,s]
    pv.power_E[h,y,s] = pv.powerMax[y,s] * ω.pv_E[h,y,s]
    # Electricity tariff
    grid.C_grid_in[h,y,s] = ω.C_grid_in[h,y,s]
    grid.C_grid_out[h,y,s] = ω.C_grid_out[h,y,s]
end
function update_investment_informations(y::Int64, s::Int64, pv::Source, liion::Liion, ω::AbstractScenarios)
    # Investment cost
    pv.C_pv[y,s] = ω.C_pv[y,s]
    liion.C_liion[y,s] = ω.C_liion[y,s]
end
# Multi-energy
function update_operation_informations(h::Int64, y::Int64, s::Int64, ld::Load, pv::Source, liion::Liion,
    h2tank::H2Tank, elyz::Electrolyzer, fc::FuelCell, tes::ThermalSto, heater::Heater, grid::Grid, ω::AbstractScenarios)
    # Energy demand and production
    ld.power_E[h,y,s] = ω.ld_E[h,y,s]
    ld.power_H[h,y,s] = ω.ld_H[h,y,s]
    pv.power_E[h,y,s] = pv.powerMax[y,s] * ω.pv_E[h,y,s]
    # Electricity tariff
    grid.C_grid_in[h,y,s] = ω.C_grid_in[h,y,s]
    grid.C_grid_out[h,y,s] = ω.C_grid_out[h,y,s]
end
function update_investment_informations(y::Int64, s::Int64, pv::Source, liion::Liion,
    h2tank::H2Tank, elyz::Electrolyzer, fc::FuelCell, tes::ThermalSto, heater::Heater, ω::AbstractScenarios)
    # Investment cost
    pv.C_pv[y,s] = ω.C_pv[y,s]
    liion.C_liion[y,s] = ω.C_liion[y,s]
    tes.C_tes[y,s] = ω.C_tes[y,s]
    h2tank.C_tank[y,s] = ω.C_tank[y,s]
    elyz.C_elyz[y,s] = ω.C_elyz[y,s]
    fc.C_fc[y,s] = ω.C_fc[y,s]
    heater.C_heater[y,s] = ω.C_heater[y,s]
end
