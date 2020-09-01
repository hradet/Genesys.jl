#
# Clustering
#_______________________________________________________________________________
# TODO: Modifier pour prendre en compte assignement
function clustering_time_block(ω::Scenarios, ntb::Int64)
    # Parameter
    ny = size(ω.ld_E,2) # numbers of years

    # interval
    tb_interval = unique([0; collect(1:ceil(Int64,ny/ntb):ny); ny])

    # time block number of years
    n_bytb = tb_interval[2:end] .- tb_interval[1:end-1]

    # Equivalent sequence of decisions by years
    σ_tb = tb_interval[2:end]

    # Time block scenarios
    ω_tb = Scenarios(ω.ld_E[:,σ_tb,:],ω.ld_H[:,σ_tb,:],ω.pv_E[:,σ_tb,:],ω.C_pv[σ_tb,:],
    ω.C_liion[σ_tb,:],ω.C_tes[σ_tb,:],ω.C_tank[σ_tb,:],ω.C_elyz[σ_tb,:],
    ω.C_fc[σ_tb,:],ω.C_heater[σ_tb,:],ω.C_grid_in[:,σ_tb,:],ω.C_grid_out[:,σ_tb,:])

    return ClusteredScenarios(ω_tb, σ_tb, n_bytb)
end
