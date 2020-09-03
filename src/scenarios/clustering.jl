mutable struct ClusteredScenarios
    clusters::NamedTuple
    σ # sequence of clusters = assignments
    nby # number of entities by cluster = counts
end

# Typical day clustering functions
function clustering_typical_day(ω::Scenarios, ntd::Int64)
    # Parameter
    nh = 24 # number of hours in one day
    nd = round(Int64, size(ω.values.ld_E,1) / nh) # number of days
    ny = size(ω.values.ld_E,2) # numbers of years
    ns = size(ω.values.ld_E,3) # numbers of scenarios

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
            ld_E = collect(reshape(ω.values.ld_E[:,y,s] / maximum(ω.values.ld_E[:,y,s]), nh, nd))
            ld_H = collect(reshape(ω.values.ld_H[:,y,s] / maximum(ω.values.ld_H[:,y,s]), nh, nd))
            pv_E = collect(reshape(ω.values.pv_E[:,y,s] / maximum(ω.values.pv_E[:,y,s]), nh, nd))

            # Combine
            data_cluster = vcat(ld_E, ld_H, pv_E)

            # Clustering
            cluster = kmeans(data_cluster, ntd)

            # Representation
            ld_E_td[:,:,y,s] = maximum(ω.values.ld_E[:,y,s]) * cluster.centers[1:nh, :]
            ld_H_td[:,:,y,s] = maximum(ω.values.ld_H[:,y,s]) * cluster.centers[nh+1:2*nh, :]
            pv_td[:,:,y,s] = maximum(ω.values.pv_E[:,y,s]) * cluster.centers[2*nh+1:3*nh, :]

            # Assignment sequence
            σ_td[:,y,s] = cluster.assignments

            # Number of assignement by cluster
            n_bytd[:,y,s] = cluster.counts
        end
    end
    # Cluster values
    clusters = (
    ld_E = ld_E_td,
    ld_H = ld_H_td,
    # Production
    pv_E = pv_td,
    # Investment costs
    C_pv = ω.values.C_pv,
    C_liion = ω.values.C_liion,
    C_tes = ω.values.C_tes,
    C_tank = ω.values.C_tank,
    C_elyz = ω.values.C_elyz,
    C_fc = ω.values.C_fc,
    C_heater = ω.values.C_heater,
    # Electricity tariff
    C_grid_in = ω.values.C_grid_in,
    C_grid_out = ω.values.C_grid_out,
    )

    return ClusteredScenarios(clusters, σ_td, n_bytd)
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

    # Cluster values
    clusters = (
    ld_E = ld_E_td,
    ld_H = ld_H_td,
    # Production
    pv_E = pv_td,
    # Investment costs
    C_pv = ω.clusters.C_pv,
    C_liion = ω.clusters.C_liion,
    C_tes = ω.clusters.C_tes,
    C_tank = ω.clusters.C_tank,
    C_elyz = ω.clusters.C_elyz,
    C_fc = ω.clusters.C_fc,
    C_heater = ω.clusters.C_heater,
    # Electricity tariff
    C_grid_in = ω.clusters.C_grid_in,
    C_grid_out = ω.clusters.C_grid_out,
    )

    return ClusteredScenarios(clusters, σ_td, n_bytd)
end
# Build scenario from a clustered scenario
function simulate_td(ω::ClusteredScenarios, y::Int64, s::Int64, horizon::Int64)

    values = (
    ld_E = reshape(ω.clusters.ld_E[:,ω.σ[:,y, s], y, s], horizon, 1),
    ld_H = reshape(ω.clusters.ld_H[:,ω.σ[:,y, s], y, s], horizon, 1),
    # Production
    pv_E = reshape(ω.clusters.pv_E[:,ω.σ[:,y, s], y, s], horizon, 1),
    # Investment costs
    C_pv = ω.clusters.C_pv,
    C_liion = ω.clusters.C_liion,
    C_tes = ω.clusters.C_tes,
    C_tank = ω.clusters.C_tank,
    C_elyz = ω.clusters.C_elyz,
    C_fc = ω.clusters.C_fc,
    C_heater = ω.clusters.C_heater,
    # Electricity tariff
    C_grid_in = ω.clusters.C_grid_in,
    C_grid_out = ω.clusters.C_grid_out,
    )

    return Scenarios(nothing, values, nothing)
end
# Time blocks function
function clustering_time_block(ω::Scenarios, ntb::Int64)
    # Parameter
    ny = size(ω.values.ld_E,2) # numbers of years

    # interval
    tb_interval = unique([0; collect(1:ceil(Int64,ny/ntb):ny); ny])

    # time block number of years
    n_bytb = tb_interval[2:end] .- tb_interval[1:end-1]

    # Equivalent sequence of decisions by years
    σ_tb = tb_interval[2:end]

    # Time block scenarios
    # Cluster values
    clusters = (
    ld_E = ω.values.ld_E[:,σ_tb,:],
    ld_H = ω.values.ld_H[:,σ_tb,:],
    # Production
    pv_E = ω.values.pv_E[:,σ_tb,:],
    # Investment costs
    C_pv = ω.values.C_pv[σ_tb,:],
    C_liion = ω.values.C_liion[σ_tb,:],
    C_tes = ω.values.C_tes[σ_tb,:],
    C_tank = ω.values.C_tank[σ_tb,:],
    C_elyz = ω.values.C_elyz[σ_tb,:],
    C_fc = ω.values.C_fc[σ_tb,:],
    C_heater = ω.values.C_heater[σ_tb,:],
    # Electricity tariff
    C_grid_in = ω.values.C_grid_in[:,σ_tb,:],
    C_grid_out = ω.values.C_grid_out[:,σ_tb,:],
    )

    return ClusteredScenarios(clusters, σ_tb, n_bytb)
end