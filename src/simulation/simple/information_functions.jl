#=
    This file includes all the funtions needed to update the operation
    and investment informations
 =#

# Operation
function update_operation_informations(h, y, s, ld, pv, liion, grid, ω)
    # Energy demand and production
    ld.power_E[h,y,s] = ω.ld_E[h,y,s]
    pv.power_E[h,y,s] = pv.powerMax[y,s] * ω.pv_E[h,y,s]
    # Electricity tariff
    grid.C_grid_in[h,y,s] = ω.C_grid_in[h,y,s]
    grid.C_grid_out[h,y,s] = ω.C_grid_out[h,y,s]
end

# Investment
function update_investment_informations(y, s, pv, liion, ω)
    # Investment cost
    pv.C_pv[y,s] = ω.C_pv[y,s]
    liion.C_liion[y,s] = ω.C_liion[y,s]
end
