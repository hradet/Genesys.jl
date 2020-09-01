#=
    This file includes all the functions needed to check the power balance
    constraints
=#

function power_balance_checking(h, y, s, ld, pv, liion, grid)
    # Electric power balance
    grid.power_E[h,y,s] = max(0. , ld.power_E[h,y,s] - pv.power_E[h,y,s] - liion.power_E[h,y,s])
end
