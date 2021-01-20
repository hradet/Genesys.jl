#=
    Scenario reduction methods
=#
abstract type AbstractScenariosReducer end
abstract type AbstractDimensionReducer end
abstract type AbstractClusteringMethod end
abstract type AbstractNormalizer end

# Manual scenario reduction
mutable struct ManualReducer <: AbstractScenariosReducer
    h::Union{UnitRange{Int64}, Int64}
    y::Union{UnitRange{Int64}, Int64}
    s::Union{UnitRange{Int64}, Int64}
    ManualReducer(; h = 1:8760, y = 1:1, s = 1:1) = new(h, y, s)
end

function reduce(reducer::ManualReducer, ω::Scenarios)
    # Parameters
    h, y, s = reducer.h, reducer.y, reducer.s
    # Demand
    ld_E = (t = ω.ld_E.t[h, y, s], power = ω.ld_E.power[h, y, s])
    ld_H = (t = ω.ld_H.t[h, y, s], power = ω.ld_H.power[h, y, s])
    # Production
    pv = (t = ω.pv.t[h, y, s], power =  ω.pv.power[h, y, s], cost =  ω.pv.cost[y, s])
    # Investment costs
    liion = (cost =  ω.liion.cost[y, s],)
    tes = (cost =  ω.tes.cost[y, s],)
    h2tank = (cost =  ω.h2tank.cost[y, s],)
    elyz = (cost =  ω.elyz.cost[y, s],)
    fc = (cost =  ω.fc.cost[y, s],)
    heater = (cost =  ω.heater.cost[y, s],)
    # Electricity tariff
    grid = (cost_in =  ω.grid.cost_in[h, y, s], cost_out =  ω.grid.cost_out[h, y, s])
    # Outputs
    ω_reduced = Scenarios(ld_E, ld_H, pv, liion, tes, h2tank, elyz, fc, heater, grid)
    probabilities = ones(length(s)) / length(s)

    return ω_reduced, probabilities
end

# Sample Average Approximation scenario reduction
mutable struct SAAReducer <: AbstractScenariosReducer
    nmontecarlo::Int64

    SAAReducer(; nmontecarlo = 100) = new(nmontecarlo)
end

function reduce(reducer::SAAReducer, ω::Scenarios{Array{DateTime,3}, Array{Float64,3}, Array{Float64,2}}; y::Int64 = 1, s::Int64 = 1)
    # Parameters
    _, ny, ns = size(ω.ld_E.power)
    # Monte carlo indices
    idx = zip(rand(y:ny, reducer.nmontecarlo), rand(s:ns, reducer.nmontecarlo))
    # Monte carlo sampling
    # Demand
    ld_E = (t = reshape(hcat([ω.ld_E.t[:, y, s] for (y,s) in idx]...),:,1,reducer.nmontecarlo), power = reshape(hcat([ω.ld_E.power[:, y, s] for (y,s) in idx]...),:,1,reducer.nmontecarlo))
    ld_H = (t = reshape(hcat([ω.ld_H.t[:, y, s] for (y,s) in idx]...),:,1,reducer.nmontecarlo), power = reshape(hcat([ω.ld_H.power[:, y, s] for (y,s) in idx]...),:,1,reducer.nmontecarlo))
    # Production
    pv = (t = reshape(hcat([ω.pv.t[:, y, s] for (y,s) in idx]...),:,1,reducer.nmontecarlo), power = reshape(hcat([ω.pv.power[:, y, s] for (y,s) in idx]...),:,1,reducer.nmontecarlo), cost = repeat(ω.pv.cost[y:y, s:s],1,reducer.nmontecarlo))
    # Electricity tariff
    grid = (cost_in = reshape(hcat([ω.grid.cost_in[:, y, s] for (y,s) in idx]...),:,1,reducer.nmontecarlo), cost_out = reshape(hcat([ω.grid.cost_out[:, y, s] for (y,s) in idx]...),:,1,reducer.nmontecarlo))
    # Investment costs
    liion = (cost =  repeat(ω.liion.cost[y:y, s:s],1,reducer.nmontecarlo),)
    tes = (cost =  repeat(ω.tes.cost[y:y, s:s],1,reducer.nmontecarlo),)
    h2tank = (cost =  repeat(ω.h2tank.cost[y:y, s:s],1,reducer.nmontecarlo),)
    elyz = (cost =  repeat(ω.elyz.cost[y:y, s:s],1,reducer.nmontecarlo),)
    fc = (cost =  repeat(ω.fc.cost[y:y, s:s],1,reducer.nmontecarlo),)
    heater = (cost =  repeat(ω.heater.cost[y:y, s:s],1,reducer.nmontecarlo),)
    # Outputs
    ω_reduced = Scenarios(ld_E, ld_H, pv, liion, tes, h2tank, elyz, fc, heater, grid)
    probabilities = ones(reducer.nmontecarlo) / reducer.nmontecarlo

    return ω_reduced, probabilities
end

# Expected value scenario reduction
mutable struct MeanValueReducer <: AbstractScenariosReducer
    MeanValueReducer() = new()
end

function reduce(reducer::MeanValueReducer, ω::Scenarios{Array{DateTime,3}, Array{Float64,3}, Array{Float64,2}}; y::Int64 = 1, s::Int64 = 1)
    # Mean value
    # Demand
    ld_E = (t = ω.ld_E.t[:, y:y, s:s], power = mean(ω.ld_E.power, dims=[2,3]))
    ld_H = (t = ω.ld_H.t[:, y:y, s:s], power = mean(ω.ld_H.power, dims=[2,3]))
    # Production
    pv = (t = ω.pv.t[:, y:y, s:s], power =  mean(ω.pv.power, dims=[2,3]), cost =  ω.pv.cost[y:y, s:s])
    # Electricity tariff
    grid = (cost_in = mean(ω.grid.cost_in, dims=[2,3]), cost_out =  mean(ω.grid.cost_out, dims=[2,3]))
    # Investment costs
    liion = (cost =  ω.liion.cost[y:y, s:s],)
    tes = (cost =  ω.tes.cost[y:y, s:s],)
    h2tank = (cost =  ω.h2tank.cost[y:y, s:s],)
    elyz = (cost =  ω.elyz.cost[y:y, s:s],)
    fc = (cost =  ω.fc.cost[y:y, s:s],)
    heater = (cost =  ω.heater.cost[y:y, s:s],)
    # Outputs
    ω_reduced = Scenarios(ld_E, ld_H, pv, liion, tes, h2tank, elyz, fc, heater, grid)
    probabilities = [1.]

    return ω_reduced, probabilities
end

# Clustering scenario reduction
mutable struct ClusteringReducer <: AbstractScenariosReducer
    normalization::AbstractNormalizer
    dimension_reduction::AbstractDimensionReducer
    clustering::AbstractClusteringMethod

    ClusteringReducer(; normalization = MinMax(),
                        dimension_reduction = UMAP(),
                        clustering = DensityBased()) = new(normalization, dimension_reduction, clustering)
end

function reduce(reducer::ClusteringReducer, ω::Scenarios{Array{DateTime,3}, Array{Float64,3}, Array{Float64,2}}; y::Int64 = 1, s::Int64 = 1)
    # Normalisation & aggregation
    norm = Genesys.normalization(reducer.normalization, ω.ld_E.power[:,2:end,:], ω.ld_H.power[:,2:end,:], ω.pv.power[:,2:end,:])#, ω.grid.cost_in[:,2:end,:], ω.grid.cost_out[:,2:end,:])
    # Dimension reduction
    embedding = Genesys.dimension_reduction(reducer.dimension_reduction, norm)
    # Clustering
    medoids, probabilities, assignments = Genesys.clustering(reducer.clustering, embedding)
    # Building reduced scenario
    # Demand
    ld_E = (t = reshape(hcat([ω.ld_E.t[:,:,s] for s in 1:ns]...)[:, medoids],:,1,reducer.clustering.n_clusters), power =  reshape(hcat([ω.ld_E.power[:,:,s] for s in 1:ns]...)[:, medoids],:,1,reducer.clustering.n_clusters))
    ld_H = (t = reshape(hcat([ω.ld_H.t[:,:,s] for s in 1:ns]...)[:, medoids],:,1,reducer.clustering.n_clusters), power = reshape(hcat([ω.ld_H.power[:,:,s] for s in 1:ns]...)[:, medoids],:,1,reducer.clustering.n_clusters))
    # Production
    pv = (t = reshape(hcat([ω.pv.t[:,:,s] for s in 1:ns]...)[:, medoids],:,1,reducer.clustering.n_clusters), power = reshape(hcat([ω.pv.power[:,:,s] for s in 1:ns]...)[:, medoids],:,1,reducer.clustering.n_clusters), cost =  repeat(ω.pv.cost[y:y, s:s],1,reducer.clustering.n_clusters))
    # Electricity tariff
    grid = (cost_in = reshape(hcat([ω.grid.cost_in[:,:,s] for s in 1:ns]...)[:, medoids],:,1,reducer.clustering.n_clusters), cost_out = reshape(hcat([ω.grid.cost_out[:,:,s] for s in 1:ns]...)[:, medoids],:,1,reducer.clustering.n_clusters))
    # Investment costs
    liion = (cost =  repeat(ω.liion.cost[y:y, s:s],1,reducer.clustering.n_clusters),)
    tes = (cost =  repeat(ω.tes.cost[y:y, s:s],1,reducer.clustering.n_clusters),)
    h2tank = (cost =  repeat(ω.h2tank.cost[y:y, s:s],1,reducer.clustering.n_clusters),)
    elyz = (cost =  repeat(ω.elyz.cost[y:y, s:s],1,reducer.clustering.n_clusters),)
    fc = (cost =  repeat(ω.fc.cost[y:y, s:s],1,reducer.clustering.n_clusters),)
    heater = (cost =  repeat(ω.heater.cost[y:y, s:s],1,reducer.clustering.n_clusters),)

    return Scenarios(ld_E, ld_H, pv, liion, tes, h2tank, elyz, fc, heater, grid), probabilities, assignments
end

# Normalization
struct MinMax <: AbstractNormalizer end

function normalization(method::MinMax, data...)
    # Dimensions
    d = size(data[1], 1)
    # Normalization
    norm = [reshape(maximum(x) == minimum(x) ? x ./ maximum(x) : (x .- minimum(x)) ./ (maximum(x) .- minimum(x)), d, :) for x in data]
    # Aggregation
    return vcat(norm...)
end

# Dimension reduction
# UMAP
mutable struct UMAP <: AbstractDimensionReducer
    n_components::Int64
    n_neighbors::Int64
    distance::Distances.SemiMetric

    UMAP(; n_components = 2, n_neighbors = 15, distance = Distances.Euclidean()) = new(n_components, n_neighbors, distance)
end

function dimension_reduction(reducer::UMAP, data::AbstractArray{Float64,2})
    # data is a d x n matrix with d dimension and n observation
    return umap(data, reducer.n_components)
end

# PCA
mutable struct PrincipalComponentAnalysis <: AbstractDimensionReducer
    n_components::Int64

    PrincipalComponentAnalysis(; n_components = 2) = new(n_components)
end

function dimension_reduction(reducer::PrincipalComponentAnalysis, data::AbstractArray{Float64,2})
    # data is a d x n matrix with d dimension and n observation
    m = MultivariateStats.fit(PCA, data, maxoutdim = reducer.n_components)
    return MultivariateStats.transform(m, data)
end

# Statistical moments
struct Moments <: AbstractDimensionReducer
    n_moments::Int64

    Moments(; n_moments = 4) = new(n_moments)
end

function dimension_reduction(reducer::Moments, data::AbstractArray{Float64,2})
    # data is a d x n matrix with d dimension and n observation
    out = hcat([moment(data[:,j], m) for j in 1:size(data,2), m in 1:reducer.n_moments])
    return permutedims(out)
end

# Clustering methods
# K-medoids
mutable struct Kmedoids <: AbstractClusteringMethod
    n_clusters::Int64
    distance::Distances.SemiMetric
    log::Bool

    Kmedoids(; n_clusters = 10, distance = Distances.Euclidean(), log = true) = new(n_clusters, distance, log)
end

function clustering(method::Kmedoids, embedding::AbstractArray{Float64,2})
    # data is a d x n matrix with d dimension and n observation
    # Distance matrix
    D = pairwise(method.distance, embedding, dims = 2)
    # Clustering
    results = kmedoids(D, method.n_clusters, display = method.log ? :iter : :none)

    return results.medoids, results.counts / sum(results.counts), results.assignments
end

# HDBSCAN
mutable struct DensityBased <: AbstractClusteringMethod
    n_clusters::Int64
    min_cluster_size::Int64

    DensityBased(; n_clusters = 10, min_cluster_size = 5) = new(n_clusters, min_cluster_size)
end

function clustering(method::DensityBased, embedding::AbstractArray{Float64,2})
    # data is a d x n matrix with d dimension and n observation
    results = HDBSCAN.hdbscan(embedding, min_cluster_size = method.min_cluster_size)
    # Retrieve vectors and of each cluster
    points = HDBSCAN.exemplars(results)
    # Find the medoids of each cluster
    dist = [pairwise(Distances.Euclidean(), permutedims(p), dims = 2 ) for p in points]
    medoids_points = [kmedoids(d, 1).medoids[1] for d in dist]
    # Get the corresponding value in the points array
    val = [p[medoids_points[i],:] for (i,p) in enumerate(points)]
    # Get the medoids from the original dataset
    medoids = [findall(v .== embedding)[1][2] for v in val]
    # Count the numbers of point in each cluster
    counts = [count(results.assignments .== c) for c in 1:maximum(results.assignments)]

    return medoids, counts / sum(counts), results.assignments
end
