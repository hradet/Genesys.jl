#=
    Scenario reduction methods
=#
abstract type AbstractScenariosReducer end
abstract type AbstractDimensionReducer end
abstract type AbstractClusteringMethod end

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
mutable struct FeatureBasedReducer <: AbstractScenariosReducer
    transformation::Union{UnionAll, Nothing}
    reduction::Union{AbstractDimensionReducer, Nothing}
    clustering::AbstractClusteringMethod

    FeatureBasedReducer(; transformation = UnitRangeTransform,
                        reduction = StatsReduction(),
                        clustering = KmedoidsClustering()) = new(transformation, reduction, clustering)
end

function reduce(reducer::FeatureBasedReducer, ω::Scenarios{Array{DateTime,3}, Array{Float64,3}, Array{Float64,2}}; y::Int64 = 1, s::Int64 = 1)
    # Parameters
    nh, ny, ns = size(ω.ld_E.power)
    # Data to be reduced : data is a tuple of d x n matrix with d dimension and n observation
    t = [reshape(ω.ld_E.t[:,2:end,:], nh, :), reshape(ω.ld_H.t[:,2:end,:], nh, :), reshape(ω.pv.t[:,2:end,:], nh, :)]
    data = [reshape(ω.ld_E.power[:,2:end,:], nh, :), reshape(ω.ld_H.power[:,2:end,:], nh, :), reshape(ω.pv.power[:,2:end,:], nh, :), reshape(ω.grid.cost_in[:,2:end,:], nh, :)]
    # Transformation
    norm = replace!.([Genesys.StatsBase.standardize(reducer.transformation, d, dims = 1) for d in data], NaN => 0.)
    # Dimension reduction
    embedding = replace!(dimension_reduction(reducer.reduction, norm), NaN => 0.)
    # Clustering
    medoids, probabilities, assignments = clustering(reducer.clustering, embedding)
    # Building reduced scenario
    # Demand
    ld_E = (t = reshape(t[1][:,medoids], nh, 1, :), power = reshape(data[1][:,medoids], nh, 1, :))
    ld_H = (t = reshape(t[2][:,medoids], nh, 1, :), power = reshape(data[2][:,medoids], nh, 1, :))
    # Production
    pv = (t = reshape(t[3][:,medoids], nh, 1, :), power = reshape(data[3][:,medoids], nh, 1, :), cost =  repeat(ω.pv.cost[y:y, s:s], 1, length(medoids)))
    # Electricity tariff
    grid = (cost_in = reshape(data[4][:,medoids], nh, 1, :), cost_out = reshape(reshape(ω.grid.cost_out[:,2:end,:], nh, :)[:,medoids], nh, 1, :))
    # Investment costs
    liion = (cost =  repeat(ω.liion.cost[y:y, s:s], 1, length(medoids)),)
    tes = (cost =  repeat(ω.tes.cost[y:y, s:s], 1, length(medoids)),)
    h2tank = (cost =  repeat(ω.h2tank.cost[y:y, s:s], 1, length(medoids)),)
    elyz = (cost =  repeat(ω.elyz.cost[y:y, s:s], 1, length(medoids)),)
    fc = (cost =  repeat(ω.fc.cost[y:y, s:s], 1, length(medoids)),)
    heater = (cost =  repeat(ω.heater.cost[y:y, s:s], 1, length(medoids)),)

    return Scenarios(ld_E, ld_H, pv, liion, tes, h2tank, elyz, fc, heater, grid), probabilities, assignments
end

# Transformation
StatsBase.standardize(DT::Nothing, X; dims=nothing, kwargs...) = X

# Dimension reduction
# UMAP
# mutable struct UMAPReduction <: AbstractDimensionReducer
#     n_components::Int64
#     n_neighbors::Int64
#     distance::Distances.SemiMetric
#
#     UMAPReduction(; n_components = 2, n_neighbors = 15, distance = Distances.Euclidean()) = new(n_components, n_neighbors, distance)
# end
#
# function dimension_reduction(reducer::UMAPReduction, data::Array{Array{Float64,2}}; aggregated::Bool=false)
#     # data is a vector of d x n matrix with d dimension and n observation
#     if aggregated
#         return umap(vcat(data...), reducer.n_components)
#     else
#         return vcat([umap(d, reducer.n_components) for d in data]...)
#     end
# end

# PCA
mutable struct PCAReduction <: AbstractDimensionReducer
    n_components::Int64

    PCAReduction(; n_components = 2) = new(n_components)
end

function dimension_reduction(reducer::PCAReduction, data::Array{Array{Float64,2}}; aggregated::Bool=false)
    # data is a vector of d x n matrix with d dimension and n observation
    if aggregated
        m = MultivariateStats.fit(PCA, vcat(data...), maxoutdim = reducer.n_components)
        return MultivariateStats.transform(m, data)
    else
        M = [MultivariateStats.fit(PCA, d, maxoutdim = reducer.n_components) for d in data]
        return vcat([MultivariateStats.transform(M[k], data[k]) for k in 1:length(data)]...)
    end
end

# Statistical moments
struct StatsReduction <: AbstractDimensionReducer end

function dimension_reduction(reducer::StatsReduction, data::Array{Array{Float64,2}}; aggregated::Bool=false)
    # data is a vector of d x n matrix with d dimension and n observation
    if aggregated
        d = vcat(data...)
        # Sum
        # s = sum(d, dims = 1)
        # Max
        # max = maximum(d, dims = 1)
        # 4 moments
        m = mean(d, dims = 1)
        v = var(d, dims = 1)
        kurt = permutedims([kurtosis(d[:,j]) for j in 1:size(d, 2)])
        skew = permutedims([skewness(d[:,j]) for j in 1:size(d, 2)])
    else
        # Sum
        # s = vcat([sum(d, dims = 1) for d in data]...)
        # Max
        # max = vcat([maximum(d, dims = 1) for d in data]...)
        # 4 moments
        m = vcat([mean(d, dims = 1) for d in data]...)
        v = vcat([var(d, dims = 1) for d in data]...)
        kurt = vcat([permutedims([kurtosis(d[:,j]) for j in 1:size(d, 2)]) for d in data]...)
        skew = vcat([permutedims([skewness(d[:,j]) for j in 1:size(d, 2)]) for d in data]...)
    end
    # Return aggregated values
    # return vcat(s, max, m, v, kurt, skew)
    return vcat(m, v, kurt, skew)
end

# No reduction
dimension_reduction(reducer::Nothing, data::Array{Array{Float64,2}}) = vcat(data...)

# Clustering methods
# K-medoids
mutable struct KmedoidsClustering <: AbstractClusteringMethod
    n_clusters::Int64
    distance::Distances.SemiMetric
    log::Bool

    KmedoidsClustering(; n_clusters = 20, distance = Distances.Euclidean(), log = true) = new(n_clusters, distance, log)
end

function clustering(method::KmedoidsClustering, embedding::AbstractArray{Float64,2})
    # data is a d x n matrix with d dimension and n observation
    # Distance matrix
    D = pairwise(method.distance, embedding, dims = 2)
    # Clustering
    results = kmedoids(D, method.n_clusters, display = method.log ? :iter : :none)

    return results.medoids, results.counts / sum(results.counts), results.assignments
end

# HDBSCAN
mutable struct HDBSCANClustering <: AbstractClusteringMethod
    n_clusters::Int64
    min_cluster_size::Int64

    HDBSCANClustering(; n_clusters = 10, min_cluster_size = 5) = new(n_clusters, min_cluster_size)
end

function clustering(method::HDBSCANClustering, embedding::AbstractArray{Float64,2})
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
