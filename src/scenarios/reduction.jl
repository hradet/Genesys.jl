#=
    Scenario reduction methods
=#
abstract type AbstractScenariosReducer end

# Manual reduction
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

# Sample Average Approximation reduction
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

# Expected value reduction
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

# Clustering reduction with kmeans
mutable struct KmeansReducer <: AbstractScenariosReducer
    ncluster::Int64

    KmeansReducer(; ncluster = 10) = new(ncluster)
end

function reduce(reducer::KmeansReducer, ω::Scenarios{Array{DateTime,3}, Array{Float64,3}, Array{Float64,2}}; y::Int64 = 1, s::Int64 = 1)
    # Parmameters
    nh = size(ω.ld_E.power,1)
    # Clustering aggregated data
    results = kmeans(ω.pv.power, ω.ld_E.power, ω.ld_H.power, ω.grid.cost_in, ω.grid.cost_out, ncluster = reducer.ncluster)
    # Demand
    ld_E = (t = repeat(ω.ld_E.t[1:nh, 1, 1], 1, 1, reducer.ncluster), power = reshape(results.centers[nh+1:2*nh,:], nh, 1, reducer.ncluster))
    ld_H = (t = repeat(ω.ld_H.t[1:nh, 1, 1], 1, 1, reducer.ncluster), power = reshape(results.centers[2*nh+1:3*nh,:], nh, 1, reducer.ncluster))
    # Production
    pv = (t = repeat(ω.pv.t[1:nh, 1, 1], 1, 1, reducer.ncluster), power = reshape(results.centers[1:nh,:], nh, 1, reducer.ncluster), cost =  repeat(ω.pv.cost[y:y, s:s],1,reducer.ncluster))
    # Electricity tariff
    grid = (cost_in = reshape(results.centers[3*nh+1:4*nh,:], nh, 1, reducer.ncluster), cost_out = reshape(results.centers[4*nh+1:5*nh,:], nh, 1, reducer.ncluster))
    # Investment costs
    liion = (cost =  repeat(ω.liion.cost[y:y, s:s],1,reducer.ncluster),)
    tes = (cost =  repeat(ω.tes.cost[y:y, s:s],1,reducer.ncluster),)
    h2tank = (cost =  repeat(ω.h2tank.cost[y:y, s:s],1,reducer.ncluster),)
    elyz = (cost =  repeat(ω.elyz.cost[y:y, s:s],1,reducer.ncluster),)
    fc = (cost =  repeat(ω.fc.cost[y:y, s:s],1,reducer.ncluster),)
    heater = (cost =  repeat(ω.heater.cost[y:y, s:s],1,reducer.ncluster),)

    return Scenarios(ld_E, ld_H, pv, liion, tes, h2tank, elyz, fc, heater, grid), results.counts / sum(results.counts),  results.assignments
end

function Clustering.kmeans(data...; ncluster::Int64 = 10)
    # Parameters
    nk = length(data)
    nh, ny, ns = size(data[1])
    # Normalization & aggregation
    data_agg = vcat([min_max_normalization(hcat([data[k][:,:,s] for s in 1:ns]...)) for k in 1:nk]...)
    # Clustering
    results = kmeans(data_agg, ncluster)
    # Denormalization
    for k in 1:nk
        results.centers[(k-1)*nh+1:k*nh,:] .= results.centers[(k-1)*nh+1:k*nh,:] .* (maximum(data[k]) .- minimum(data[k])) .+ minimum(data[k])
    end
    return results
end

# Clustering reduction with kmedoids
mutable struct KmedoidsReducer <: AbstractScenariosReducer
    ncluster::Int64
    distance

    KmedoidsReducer(; ncluster = 10, distance = Distances.Euclidean()) = new(ncluster, distance)
end

function reduce(reducer::KmedoidsReducer, ω::Scenarios{Array{DateTime,3}, Array{Float64,3}, Array{Float64,2}}; y::Int64 = 1, s::Int64 = 1)
    # Parmameters
    ns = size(ω.ld_E.power,3)
    # Clustering aggregated data
    results = kmedoids(ω.pv.power, ω.ld_E.power, ω.ld_H.power, ω.grid.cost_in, ω.grid.cost_out, ncluster = reducer.ncluster, distance = reducer.distance, display = :final)
    # Demand
    ld_E = (t = reshape(hcat([ω.ld_E.t[:,:,s] for s in 1:ns]...)[:, results.medoids],:,1,reducer.ncluster), power =  reshape(hcat([ω.ld_E.power[:,:,s] for s in 1:ns]...)[:, results.medoids],:,1,reducer.ncluster))
    ld_H = (t = reshape(hcat([ω.ld_H.t[:,:,s] for s in 1:ns]...)[:, results.medoids],:,1,reducer.ncluster), power = reshape(hcat([ω.ld_H.power[:,:,s] for s in 1:ns]...)[:, results.medoids],:,1,reducer.ncluster))
    # Production
    pv = (t = reshape(hcat([ω.pv.t[:,:,s] for s in 1:ns]...)[:, results.medoids],:,1,reducer.ncluster), power = reshape(hcat([ω.pv.power[:,:,s] for s in 1:ns]...)[:, results.medoids],:,1,reducer.ncluster), cost =  repeat(ω.pv.cost[y:y, s:s],1,reducer.ncluster))
    # Electricity tariff
    grid = (cost_in = reshape(hcat([ω.grid.cost_in[:,:,s] for s in 1:ns]...)[:, results.medoids],:,1,reducer.ncluster), cost_out = reshape(hcat([ω.grid.cost_out[:,:,s] for s in 1:ns]...)[:, results.medoids],:,1,reducer.ncluster))
    # Investment costs
    liion = (cost =  repeat(ω.liion.cost[y:y, s:s],1,reducer.ncluster),)
    tes = (cost =  repeat(ω.tes.cost[y:y, s:s],1,reducer.ncluster),)
    h2tank = (cost =  repeat(ω.h2tank.cost[y:y, s:s],1,reducer.ncluster),)
    elyz = (cost =  repeat(ω.elyz.cost[y:y, s:s],1,reducer.ncluster),)
    fc = (cost =  repeat(ω.fc.cost[y:y, s:s],1,reducer.ncluster),)
    heater = (cost =  repeat(ω.heater.cost[y:y, s:s],1,reducer.ncluster),)

    return Scenarios(ld_E, ld_H, pv, liion, tes, h2tank, elyz, fc, heater, grid), results.counts / sum(results.counts), results.assignments
end

function Clustering.kmedoids(data...; ncluster::Int64 = 10, distance = Distances.Euclidean(), display::Symbol = :iter)
    # Parameters
    nk = length(data)
    nh, ny, ns = size(data[1])
    # Normalization & aggregation
    data_agg = vcat([min_max_normalization(hcat([data[k][:,:,s] for s in 1:ns]...)) for k in 1:nk]...)
    # Compute euclidean distance matrix
    dist = pairwise(distance, data_agg, dims = 2 )
    # Clustering
    results = kmedoids(dist, ncluster, display=display)
    return results
end

# Clustering reduction with kmedoids
mutable struct FeatureBasedReducer <: AbstractScenariosReducer
    ncluster::Int64
    distance

    FeatureBasedReducer(; ncluster = 10, distance = Distances.Euclidean()) = new(ncluster, distance)
end

function reduce(reducer::FeatureBasedReducer, ω::Scenarios{Array{DateTime,3}, Array{Float64,3}, Array{Float64,2}}; y::Int64 = 1, s::Int64 = 1)
    # Parmameters
    ny, ns = size(ω.ld_E.power,2), size(ω.ld_E.power,3)
    # Features extraction - mean, variance, skewness, kurtosis
    pv_features = vcat(mean(ω.pv.power, dims=1), var(ω.pv.power, dims=1), reshape([skewness(ω.pv.power[:,y,s]) for y in 1:ny, s in 1:ns],1,ny,ns), reshape([kurtosis(ω.pv.power[:,y,s]) for y in 1:ny, s in 1:ns],1,ny,ns))
    ld_E_features = vcat(mean(ω.ld_E.power, dims=1), var(ω.ld_E.power, dims=1), reshape([skewness(ω.ld_E.power[:,y,s]) for y in 1:ny, s in 1:ns],1,ny,ns), reshape([kurtosis(ω.ld_E.power[:,y,s]) for y in 1:ny, s in 1:ns],1,ny,ns))
    ld_H_features = vcat(mean(ω.ld_H.power, dims=1), var(ω.ld_H.power, dims=1), reshape([skewness(ω.ld_H.power[:,y,s]) for y in 1:ny, s in 1:ns],1,ny,ns), reshape([kurtosis(ω.ld_H.power[:,y,s]) for y in 1:ny, s in 1:ns],1,ny,ns))
    cost_in_features = vcat(mean(ω.grid.cost_in, dims=1), var(ω.grid.cost_in, dims=1), reshape([skewness(ω.grid.cost_in[:,y,s]) for y in 1:ny, s in 1:ns],1,ny,ns), reshape([kurtosis(ω.grid.cost_in[:,y,s]) for y in 1:ny, s in 1:ns],1,ny,ns))
    cost_out_features = vcat(mean(ω.grid.cost_out, dims=1), var(ω.grid.cost_out, dims=1), reshape([skewness(ω.grid.cost_out[:,y,s]) for y in 1:ny, s in 1:ns],1,ny,ns), reshape([kurtosis(ω.grid.cost_out[:,y,s]) for y in 1:ny, s in 1:ns],1,ny,ns))
    # Clustering aggregated data
    results = kmedoids(pv_features, ld_E_features, ld_H_features, cost_in_features, cost_out_features, ncluster = reducer.ncluster, distance = reducer.distance, display = :final)
    # Demand
    ld_E = (t = reshape(hcat([ω.ld_E.t[:,:,s] for s in 1:ns]...)[:, results.medoids],:,1,reducer.ncluster), power =  reshape(hcat([ω.ld_E.power[:,:,s] for s in 1:ns]...)[:, results.medoids],:,1,reducer.ncluster))
    ld_H = (t = reshape(hcat([ω.ld_H.t[:,:,s] for s in 1:ns]...)[:, results.medoids],:,1,reducer.ncluster), power = reshape(hcat([ω.ld_H.power[:,:,s] for s in 1:ns]...)[:, results.medoids],:,1,reducer.ncluster))
    # Production
    pv = (t = reshape(hcat([ω.pv.t[:,:,s] for s in 1:ns]...)[:, results.medoids],:,1,reducer.ncluster), power = reshape(hcat([ω.pv.power[:,:,s] for s in 1:ns]...)[:, results.medoids],:,1,reducer.ncluster), cost =  repeat(ω.pv.cost[y:y, s:s],1,reducer.ncluster))
    # Electricity tariff
    grid = (cost_in = reshape(hcat([ω.grid.cost_in[:,:,s] for s in 1:ns]...)[:, results.medoids],:,1,reducer.ncluster), cost_out = reshape(hcat([ω.grid.cost_out[:,:,s] for s in 1:ns]...)[:, results.medoids],:,1,reducer.ncluster))
    # Investment costs
    liion = (cost =  repeat(ω.liion.cost[y:y, s:s],1,reducer.ncluster),)
    tes = (cost =  repeat(ω.tes.cost[y:y, s:s],1,reducer.ncluster),)
    h2tank = (cost =  repeat(ω.h2tank.cost[y:y, s:s],1,reducer.ncluster),)
    elyz = (cost =  repeat(ω.elyz.cost[y:y, s:s],1,reducer.ncluster),)
    fc = (cost =  repeat(ω.fc.cost[y:y, s:s],1,reducer.ncluster),)
    heater = (cost =  repeat(ω.heater.cost[y:y, s:s],1,reducer.ncluster),)

    return Scenarios(ld_E, ld_H, pv, liion, tes, h2tank, elyz, fc, heater, grid), results.counts / sum(results.counts), results.assignments
end
