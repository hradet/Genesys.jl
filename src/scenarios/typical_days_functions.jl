#
# Clustering
#_______________________________________________________________________________

function clustering_typical_day(ω::Scenarios, ntd::Int64)
    # Parameter
    nh = 24 # number of hours in one day
    nd = round(Int64, size(ω.ld_E,1) / nh) # number of days
    ny = size(ω.ld_E,2) # numbers of years
    ns = size(ω.ld_E,3) # numbers of scenarios

    # Pre allocate
    ld_E_td = zeros(nh, ntd, ny, ns)
    ld_H_td = zeros(nh, ntd, ny, ns)
    pv_td = zeros(nh, ntd, ny, ns)
    σ_td = zeros(Int64, nd, ny, ns)
    n_bytd = zeros(Int64, ntd, ny, ns)

    # For each year of each scenario
    for s in 1:ns
        for y in 1:ny
            # Normalization
            ld_E = collect(reshape(ω.ld_E[:,y,s] / maximum(ω.ld_E[:,y,s]), nh, nd))
            ld_H = collect(reshape(ω.ld_H[:,y,s] / maximum(ω.ld_H[:,y,s]), nh, nd))
            pv_E = collect(reshape(ω.pv_E[:,y,s] / maximum(ω.pv_E[:,y,s]), nh, nd))

            # Combine
            data_cluster = vcat(ld_E, ld_H, pv_E)

            # Clustering
            cluster = kmeans(data_cluster, ntd)

            # Representation
            ld_E_td[:,:,y,s] = maximum(ω.ld_E[:,y,s]) * cluster.centers[1:nh, :]
            ld_H_td[:,:,y,s] = maximum(ω.ld_H[:,y,s]) * cluster.centers[nh+1:2*nh, :]
            pv_td[:,:,y,s] = maximum(ω.pv_E[:,y,s]) * cluster.centers[2*nh+1:3*nh, :]

            # Assignment sequence
            σ_td[:,y,s] = cluster.assignments

            # Number of assignement by cluster
            n_bytd[:,y,s] = cluster.counts
        end
    end
    # Tupical day scenarios
    ω_td = Scenarios(ld_E_td, ld_H_td, pv_td,ω.C_pv,ω.C_liion,ω.C_tes,ω.C_tank,ω.C_elyz,ω.C_fc,ω.C_heater,ω.C_grid_in,ω.C_grid_out)

    return ClusteredScenarios(ω_td, σ_td, n_bytd)
end
function clustering_typical_day(ω::ClusteredScenarios, ntd::Int64)
    # Parameter
    nh = 24 # number of hours in one day
    nd = round(Int64, size(ω.clusters.ld_E,1) / nh) # number of days
    ny = size(ω.clusters.ld_E,2) # numbers of years
    ns = size(ω.clusters.ld_E,3) # numbers of scenarios

    # Pre allocate
    ld_E_td = zeros(nh, ntd, ny, ns)
    ld_H_td = zeros(nh, ntd, ny, ns)
    pv_td = zeros(nh, ntd, ny, ns)
    σ_td = zeros(Int64, nd, ny, ns)
    n_bytd = zeros(Int64, ntd, ny, ns)

    # For each year of each scenario
    for s in 1:ns
        for y in 1:ny
            # Normalization
            ld_E = collect(reshape(ω.clusters.ld_E[:,y,s] / maximum(ω.clusters.ld_E[:,y,s]), nh, nd))
            ld_H = collect(reshape(ω.clusters.ld_H[:,y,s] / maximum(ω.clusters.ld_H[:,y,s]), nh, nd))
            pv_E = collect(reshape(ω.clusters.pv_E[:,y,s] / maximum(ω.clusters.pv_E[:,y,s]), nh, nd))

            # Combine
            data_cluster = vcat(ld_E, ld_H, pv_E)

            # Clustering
            cluster = kmeans(data_cluster, ntd)

            # Representation
            ld_E_td[:,:,y,s] = maximum(ω.clusters.ld_E[:,y,s]) * cluster.centers[1:nh, :]
            ld_H_td[:,:,y,s] = maximum(ω.clusters.ld_H[:,y,s]) * cluster.centers[nh+1:2*nh, :]
            pv_td[:,:,y,s] = maximum(ω.clusters.pv_E[:,y,s]) * cluster.centers[2*nh+1:3*nh, :]

            # Assignment sequence
            σ_td[:,y,s] = cluster.assignments

            # Number of assignement by cluster
            n_bytd[:,y,s] = cluster.counts
        end
    end

    # Tupical day scenarios
    ω_td = Scenarios(ld_E_td, ld_H_td, pv_td,ω.clusters.C_pv,ω.clusters.C_liion,
    ω.clusters.C_tes,ω.clusters.C_tank,ω.clusters.C_elyz,ω.clusters.C_fc,
    ω.clusters.C_heater,ω.clusters.C_grid_in,ω.clusters.C_grid_out)

    return ClusteredScenarios(ω_td, σ_td, n_bytd)
end

#
# Simulate
#_______________________________________________________________________________
function simulate_td(ω::ClusteredScenarios, y::Int64, s::Int64, horizon::Int64)

    ld_E_td = reshape(ω.clusters.ld_E[:,ω.σ[:,y, s], y, s], horizon, 1)
    ld_H_td = reshape(ω.clusters.ld_H[:,ω.σ[:,y, s], y, s], horizon, 1)
    pv_td = reshape(ω.clusters.pv_E[:,ω.σ[:,y, s], y, s], horizon, 1)

    return Scenarios(ld_E_td, ld_H_td, pv_td,ω.clusters.C_pv,ω.clusters.C_liion,
    ω.clusters.C_tes,ω.clusters.C_tank,ω.clusters.C_elyz,ω.clusters.C_fc,
    ω.clusters.C_heater,ω.clusters.C_grid_in,ω.clusters.C_grid_out)
end
