mutable struct ClusteredScenarios
    ω::Scenarios
    assignements::Vector{Int64}
    counts::Vector{Int64}
end


function clustering(data...; ncluster::Int64 = 1, algo = "kmeans")
    # Parameters
    nk = length(data)
    nh, ny, ns = size(data[1])
    # Normalization
    data_n = data ./ maximum.(data)
    # Aggregation
    data_agg = vcat([hcat([data_n[k][:,:,s] for s in 1:ns]...) for k in 1:nk]...)
    if algo == "kmeans"
        # Clustering
        results = kmeans(data_agg, ncluster)
        # Denormalization
        for k in 1:nk
            results.centers[(k-1)*nh+1:k*nh,:] .*= maximum(data[k])
        end
    elseif algo == "kmedoids"
        # Compute euclidean distance matrix
        dist = pairwise(Euclidean(), data_agg, dims = 2 )
        # Clustering
        results = kmedoids(dist, ncluster)
    end

    return results
end

# Typical day clustering functions
function clustering_typical_day(ω::Scenarios, ntd::Int64)
    # Parameter
    nh = 24 # number of hours in one day
    nd = round(Int64, size(ω.ld_E.power,1) / nh) # number of days
    ny = size(ω.ld_E.power,2) # numbers of years
    ns = size(ω.ld_E.power,3) # numbers of scenarios

    # Pre allocate
    ld_E_td = zeros(nh, ntd, ny, ns)
    ld_H_td = zeros(nh, ntd, ny, ns)
    pv_td = zeros(nh, ntd, ny, ns)
    assignements = zeros(Int64, nd, ny, ns)
    n_bytd = zeros(Int64, ntd, ny, ns)

    # For each year of each scenario
    for s in 1:ns
        for y in 1:ny
            # Normalization
            ld_E = collect(reshape(ω.ld_E.power[:,y,s] / maximum(ω.ld_E.power[:,y,s]), nh, nd))
            ld_H = collect(reshape(ω.ld_H.power[:,y,s] / maximum(ω.ld_H.power[:,y,s]), nh, nd))
            pv = collect(reshape(ω.pv.power[:,y,s] / maximum(ω.pv.power[:,y,s]), nh, nd))

            # Combine
            data_cluster = vcat(ld_E, ld_H, pv)

            # Clustering
            cluster = kmeans(data_cluster, ntd)

            # Representation
            ld_E_td[:,:,y,s] = maximum(ω.ld_E.power[:,y,s]) * cluster.centers[1:nh, :]
            ld_H_td[:,:,y,s] = maximum(ω.ld_H.power[:,y,s]) * cluster.centers[nh+1:2*nh, :]
            pv_td[:,:,y,s] = maximum(ω.pv.power[:,y,s]) * cluster.centers[2*nh+1:3*nh, :]

            # Assignment sequence
            assignements[:,y,s] = cluster.assignments

            # Number of assignement by cluster
            counts[:,y,s] = cluster.counts
        end
    end
    # Cluster values
    clusters = (
    ld_E = ld_E_td,
    ld_H = ld_H_td,
    # Production
    pv = pv_td,
    # Investment costs
    C_pv = ω.C_pv,
    C_liion = ω.C_liion,
    C_tes = ω.C_tes,
    C_tank = ω.C_tank,
    C_elyz = ω.C_elyz,
    C_fc = ω.C_fc,
    C_heater = ω.C_heater,
    # Electricity tariff
    C_grid_in = ω.C_grid_in,
    C_grid_out = ω.C_grid_out,
    )

    return ClusteredScenarios(clusters, assignements, counts)
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
            pv_E = collect(reshape(ω.clusters.pv[:,y,s] / maximum(ω.clusters.pv[:,y,s]), nh, nd))

            # Combine
            data_cluster = vcat(ld_E, ld_H, pv)

            # Clustering
            cluster = kmeans(data_cluster, ntd)

            # Representation
            ld_E_td[:,:,y,s] = maximum(ω.clusters.ld_E[:,y,s]) * cluster.centers[1:nh, :]
            ld_H_td[:,:,y,s] = maximum(ω.clusters.ld_H[:,y,s]) * cluster.centers[nh+1:2*nh, :]
            pv_td[:,:,y,s] = maximum(ω.clusters.pv[:,y,s]) * cluster.centers[2*nh+1:3*nh, :]

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
    pv = pv_td,
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
    pv = reshape(ω.clusters.pv[:,ω.σ[:,y, s], y, s], horizon, 1),
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
    ny = size(ω.ld_E.power,2) # numbers of years

    # interval
    tb_interval = unique([0; collect(1:ceil(Int64,ny/ntb):ny); ny])

    # time block number of years
    n_bytb = tb_interval[2:end] .- tb_interval[1:end-1]

    # Equivalent sequence of decisions by years
    σ_tb = tb_interval[2:end]

    # Time block scenarios
    # Cluster values
    clusters = (
    ld_E = ω.ld_E.power[:,σ_tb,:],
    ld_H = ω.ld_H.power[:,σ_tb,:],
    # Production
    pv = ω.pv.power[:,σ_tb,:],
    # Investment costs
    C_pv = ω.C_pv[σ_tb,:],
    C_liion = ω.C_liion[σ_tb,:],
    C_tes = ω.C_tes[σ_tb,:],
    C_tank = ω.C_tank[σ_tb,:],
    C_elyz = ω.C_elyz[σ_tb,:],
    C_fc = ω.C_fc[σ_tb,:],
    C_heater = ω.C_heater[σ_tb,:],
    # Electricity tariff
    C_grid_in = ω.C_grid_in[:,σ_tb,:],
    C_grid_out = ω.C_grid_out[:,σ_tb,:],
    )

    return ClusteredScenarios(clusters, σ_tb, n_bytb)
end
