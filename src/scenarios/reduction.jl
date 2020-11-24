#=
    Scenario reduction methods
=#
#TODO différencier multiyear scenarios de ceux d'une année ??
abstract type AbstractScenarioReduction end

# Manual reduction
mutable struct ManualReduction <: AbstractScenarioReduction
    ManualReduction() = new()
end

function reduce(reducer::ManualReduction, ω::Scenarios; h::Union{UnitRange{Int64}, Int64}= 1:8760, y::Union{UnitRange{Int64}, Int64}= 1, s::Union{UnitRange{Int64}, Int64}= 1)
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
    probabilities = ones(length(y), length(s)) / length(y)

    return ω_reduced, probabilities
end

# Sample Average Approximation reduction
mutable struct SAAReduction <: AbstractScenarioReduction
    nmontecarlo::Int64

    SAAReduction(; nmontecarlo = 100) = new(nmontecarlo)
end

function reduce(reducer::SAAReduction, ω::Scenarios; y::Int64 = 1, s::Int64 = 1)
    # Parameters
    _, ny, ns = size(ω.ld_E.power)
    # Monte carlo indices
    idx = zip(rand(y:ny, reducer.nmontecarlo), rand(s:ns, reducer.nmontecarlo))
    # Monte carlo sampling
    # Demand
    ld_E = (t = hcat([ω.ld_E.t[:, y, s] for (y,s) in idx]...), power = hcat([ω.ld_E.power[:, y, s] for (y,s) in idx]...))
    ld_H = (t = hcat([ω.ld_H.t[:, y, s] for (y,s) in idx]...), power = hcat([ω.ld_H.power[:, y, s] for (y,s) in idx]...))
    # Production
    pv = (t = hcat([ω.pv.t[:, y, s] for (y,s) in idx]...), power = hcat([ω.pv.power[:, y, s] for (y,s) in idx]...), cost = ω.pv.cost[y, s])
    # Electricity tariff
    grid = (cost_in = hcat([ω.grid.cost_in[:, y, s] for (y,s) in idx]...), cost_out = hcat([ω.grid.cost_out[:, y, s] for (y,s) in idx]...))
    # Investment costs
    liion = (cost =  ω.liion.cost[y, s],)
    tes = (cost =  ω.tes.cost[y, s],)
    h2tank = (cost =  ω.h2tank.cost[y, s],)
    elyz = (cost =  ω.elyz.cost[y, s],)
    fc = (cost =  ω.fc.cost[y, s],)
    heater = (cost =  ω.heater.cost[y, s],)
    # Outputs
    ω_reduced = Scenarios(ld_E, ld_H, pv, liion, tes, h2tank, elyz, fc, heater, grid)
    probabilities = ones(reducer.nmontecarlo) / reducer.nmontecarlo

    return ω_reduced, probabilities
end

# Expected value reduction
mutable struct ExpectedValueReduction <: AbstractScenarioReduction
    ExpectedValueReduction() = new()
end

function reduce(reducer::ExpectedValueReduction, ω::Scenarios; y::Int64 = 1, s::Int64 = 1)
    # Mean value
    # Demand
    ld_E = (t = ω.ld_E.t[:, 1, 1], power = mean(ω.ld_E.power[:, :, :], dims=[2,3]))
    ld_H = (t = ω.ld_H.t[:, 1, 1], power = mean(ω.ld_H.power[:, :, :], dims=[2,3]))
    # Production
    pv = (t = ω.pv.t[:, 1, 1], power =  mean(ω.pv.power[:, :, :], dims=[2,3]), cost =  ω.pv.cost[y, s])
    # Electricity tariff
    grid = (cost_in = mean(ω.grid.cost_in[:, :, :], dims=[2,3]), cost_out =  mean(ω.grid.cost_out[:, :, :], dims=[2,3]))
    # Investment costs
    liion = (cost =  ω.liion.cost[y, s],)
    tes = (cost =  ω.tes.cost[y, s],)
    h2tank = (cost =  ω.h2tank.cost[y, s],)
    elyz = (cost =  ω.elyz.cost[y, s],)
    fc = (cost =  ω.fc.cost[y, s],)
    heater = (cost =  ω.heater.cost[y, s],)
    # Outputs
    ω_reduced = Scenarios(ld_E, ld_H, pv, liion, tes, h2tank, elyz, fc, heater, grid)
    probabilities = [1.]

    return ω_reduced, probabilities
end

# Clustering reduction with kmeans
mutable struct KmeansReduction <: AbstractScenarioReduction
    ncluster::Int64
    typical_days::Bool

    KmeansReduction(; ncluster = 10, typical_days = false) = new(ncluster, typical_days)
end

function reduce(reducer::KmeansReduction, ω::Scenarios; y::Int64 = 1, s::Int64 = 1)
    if reducer.typical_days
        # Parmameters
        nh = 24
        # Clustering aggregated data
        results = kmeans(reshape(ω.pv.power[:,y,s], nh, :, 1), reshape(ω.ld_E.power[:,y,s], nh, :, 1), reshape(ω.ld_H.power[:,y,s], nh, :, 1), reshape(ω.grid.cost_in[:,y,s], nh, :, 1), reshape(ω.grid.cost_out[:,y,s], nh, :, 1), ncluster = reducer.ncluster)
    else
        # Parmameters
        nh = size(ω.ld_E.power,1)
        # Clustering aggregated data
        results = kmeans(ω.pv.power, ω.ld_E.power, ω.ld_H.power, ω.grid.cost_in, ω.grid.cost_out, ncluster = reducer.ncluster)
    end
    # Demand
    ld_E = (t = repeat(ω.ld_E.t[1:nh, y, s], 1, reducer.ncluster), power = results.centers[nh+1:2*nh,:])
    ld_H = (t = repeat(ω.ld_H.t[1:nh, y, s], 1, reducer.ncluster), power = results.centers[2*nh+1:3*nh,:])
    # Production
    pv = (t = repeat(ω.pv.t[1:nh, y, s], 1, reducer.ncluster), power = results.centers[1:nh,:], cost =  ω.pv.cost[y, s])
    # Electricity tariff
    grid = (cost_in = results.centers[3*nh+1:4*nh,:], cost_out = results.centers[4*nh+1:5*nh,:])
    # Investment costs
    liion = (cost =  ω.liion.cost[y, s],)
    tes = (cost =  ω.tes.cost[y, s],)
    h2tank = (cost =  ω.h2tank.cost[y, s],)
    elyz = (cost =  ω.elyz.cost[y, s],)
    fc = (cost =  ω.fc.cost[y, s],)
    heater = (cost =  ω.heater.cost[y, s],)

    return Scenarios(ld_E, ld_H, pv, liion, tes, h2tank, elyz, fc, heater, grid), results.counts, results.assignments
end

function Clustering.kmeans(data...; ncluster::Int64 = 10)
    # Parameters
    nk = length(data)
    nh, ny, ns = size(data[1])
    # Normalization
    data_n = data ./ maximum.(data)
    # Aggregation
    data_agg = vcat([hcat([data_n[k][:,:,s] for s in 1:ns]...) for k in 1:nk]...)
    # Clustering
    results = kmeans(data_agg, ncluster)
    # Denormalization
    for k in 1:nk
        results.centers[(k-1)*nh+1:k*nh,:] .*= maximum(data[k])
    end
    return results
end

# Clustering reduction with kmedoids
mutable struct KmedoidsReduction <: AbstractScenarioReduction
    ncluster::Int64
    distance
    typical_days::Bool

    KmedoidsReduction(; ncluster = 10, distance = Distances.Euclidean(), typical_days = false) = new(ncluster, distance, typical_days)
end

function reduce(reducer::KmedoidsReduction, ω::Scenarios; y::Int64 = 1, s::Int64 = 1)
    if reducer.typical_days
        # Parmameters
        nh, ns = 24, 1
        # Clustering aggregated data
        results = kmedoids(reshape(ω.pv.power[:,y,s], nh, :, ns), reshape(ω.ld_E.power[:,y,s], nh, :, ns), reshape(ω.ld_H.power[:,y,s], nh, :, ns), reshape(ω.grid.cost_in[:,y,s], nh, :, ns), reshape(ω.grid.cost_out[:,y,s], nh, :, ns), ncluster = reducer.ncluster, distance = reducer.distance)
        # Demand
        ld_E = (t = hcat([reshape(ω.ld_E.t[:,y,s], nh, :, ns)[:,:,s] for s in 1:ns]...)[:, results.medoids], power =  hcat([reshape(ω.ld_E.power[:,y,s], nh, :, ns)[:,:,s] for s in 1:ns]...)[:, results.medoids])
        ld_H = (t = hcat([reshape(ω.ld_H.t[:,y,s], nh, :, ns)[:,:,s] for s in 1:ns]...)[:, results.medoids], power = hcat([reshape(ω.ld_H.power[:,y,s], nh, :, ns)[:,:,s] for s in 1:ns]...)[:, results.medoids])
        # Production
        pv = (t = hcat([reshape(ω.pv.t[:,y,s], nh, :, ns)[:,:,s] for s in 1:ns]...)[:, results.medoids], power = hcat([reshape(ω.pv.power[:,y,s], nh, :, ns)[:,:,s] for s in 1:ns]...)[:, results.medoids], cost =  ω.pv.cost[y, s])
        # Electricity tariff
        grid = (cost_in = hcat([reshape(ω.grid.cost_in[:,y,s], nh, :, ns)[:,:,s] for s in 1:ns]...)[:, results.medoids], cost_out = hcat([reshape(ω.grid.cost_out[:,y,s], nh, :, ns)[:,:,s] for s in 1:ns]...)[:, results.medoids])
    else
        # Parmameters
        ns = size(ω.ld_E.power, 3)
        # Clustering aggregated data
        results = kmedoids(ω.pv.power, ω.ld_E.power, ω.ld_H.power, ω.grid.cost_in, ω.grid.cost_out, ncluster = reducer.ncluster, distance = reducer.distance)
        # Demand
        ld_E = (t = hcat([ω.ld_E.t[:,:,s] for s in 1:ns]...)[:, results.medoids], power =  hcat([ω.ld_E.power[:,:,s] for s in 1:ns]...)[:, results.medoids])
        ld_H = (t = hcat([ω.ld_H.t[:,:,s] for s in 1:ns]...)[:, results.medoids], power = hcat([ω.ld_H.power[:,:,s] for s in 1:ns]...)[:, results.medoids])
        # Production
        pv = (t = hcat([ω.pv.t[:,:,s] for s in 1:ns]...)[:, results.medoids], power = hcat([ω.pv.power[:,:,s] for s in 1:ns]...)[:, results.medoids], cost =  ω.pv.cost[y, s])
        # Electricity tariff
        grid = (cost_in = hcat([ω.grid.cost_in[:,:,s] for s in 1:ns]...)[:, results.medoids], cost_out = hcat([ω.grid.cost_out[:,:,s] for s in 1:ns]...)[:, results.medoids])
    end
    # Investment costs
    liion = (cost =  ω.liion.cost[y, s],)
    tes = (cost =  ω.tes.cost[y, s],)
    h2tank = (cost =  ω.h2tank.cost[y, s],)
    elyz = (cost =  ω.elyz.cost[y, s],)
    fc = (cost =  ω.fc.cost[y, s],)
    heater = (cost =  ω.heater.cost[y, s],)

    return Scenarios(ld_E, ld_H, pv, liion, tes, h2tank, elyz, fc, heater, grid), results.counts, results.assignments
end

function Clustering.kmedoids(data...; ncluster::Int64 = 10, distance = Distances.Euclidean())
    # Parameters
    nk = length(data)
    nh, ny, ns = size(data[1])
    # Normalization
    data_n = data ./ maximum.(data)
    # Aggregation
    data_agg = vcat([hcat([data_n[k][:,:,s] for s in 1:ns]...) for k in 1:nk]...)
    # Compute euclidean distance matrix
    dist = pairwise(distance, data_agg, dims = 2 )
    # Clustering
    results = kmedoids(dist, ncluster)
    return results
end
