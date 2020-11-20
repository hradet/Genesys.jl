#=
    Scenario reduction functions
=#
mutable struct Scenarios{A,B,C,D} <: AbstractScenarios
    # Demand
    ld_E::A
    ld_H::A
    # Production
    pv::B
    # Investment costs
    liion::C
    tes::C
    h2tank::C
    elyz::C
    fc::C
    heater::C
    # Electricity tariff
    grid::D
end

# Constructor
function Scenarios(d::Dict{}, h::Union{UnitRange{Int64}, Int64}, y::Union{UnitRange{Int64}, Int64}, s::Union{UnitRange{Int64}, Int64})
    # Demand
    ld_E = (t = d["ld_E"].t[h, y, s], power = d["ld_E"].power[h, y, s])
    ld_H = (t = d["ld_H"].t[h, y, s], power = d["ld_H"].power[h, y, s])
    # Production
    pv = (t = d["pv"].t[h, y, s], power = d["pv"].power[h, y, s], cost = d["pv"].cost[y, s])
    # Investment costs
    liion = (cost = d["liion"].cost[y, s],)
    tes = (cost = d["tes"].cost[y, s],)
    h2tank = (cost = d["h2tank"].cost[y, s],)
    elyz = (cost = d["elyz"].cost[y, s],)
    fc = (cost = d["fc"].cost[y, s],)
    heater = (cost = d["heater"].cost[y, s],)
    # Electricity tariff
    grid = (cost_in = d["grid"].cost_in[h, y, s], cost_out = d["grid"].cost_out[h, y, s])

    return Scenarios(ld_E, ld_H, pv, liion, tes, h2tank, elyz, fc, heater, grid)
end

# Scenario Reduction functions
function scenarios_reduction(ω::Scenarios, h::Union{UnitRange{Int64}, Int64}, y::Union{UnitRange{Int64}, Int64}, s::Union{UnitRange{Int64}, Int64})
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

    return Scenarios(ld_E, ld_H, pv, liion, tes, h2tank, elyz, fc, heater, grid)
end

function scenarios_reduction_SAA(ω::Scenarios, h::Union{UnitRange{Int64}, Int64}, y::Int64, s::Int64, Nmc::Int64)
    ny = size(ω.ld_E.power, 2)
    ns = size(ω.ld_E.power, 3)

    # Monte carlo indices
    idx = zip(rand(y:ny, Nmc), rand(1:ns, Nmc))

    # Monte carlo sampling
    # Demand
    ld_E = (t = hcat([ω.ld_E.t[h, y, s] for (y,s) in idx]...), power = hcat([ω.ld_E.power[h, y, s] for (y,s) in idx]...))
    ld_H = (t = hcat([ω.ld_H.t[h, y, s] for (y,s) in idx]...), power = hcat([ω.ld_H.power[h, y, s] for (y,s) in idx]...))
    # Production
    pv = (t = hcat([ω.pv.t[h, y, s] for (y,s) in idx]...), power = hcat([ω.pv.power[h, y, s] for (y,s) in idx]...), cost = ω.pv.cost[y, s])
    # Electricity tariff
    grid = (cost_in = hcat([ω.grid.cost_in[h, y, s] for (y,s) in idx]...), cost_out = hcat([ω.grid.cost_out[h, y, s] for (y,s) in idx]...))

    # Investment costs
    liion = (cost =  ω.liion.cost[y, s],)
    tes = (cost =  ω.tes.cost[y, s],)
    h2tank = (cost =  ω.h2tank.cost[y, s],)
    elyz = (cost =  ω.elyz.cost[y, s],)
    fc = (cost =  ω.fc.cost[y, s],)
    heater = (cost =  ω.heater.cost[y, s],)

    return Scenarios(ld_E, ld_H, pv, liion, tes, h2tank, elyz, fc, heater, grid)
end

function scenarios_reduction_mean(ω::Scenarios, h::Union{UnitRange{Int64}, Int64}, y::Int64, s::Int64)
    # Mean value over y and s
    # Demand
    ld_E = (t = ω.ld_E.t[h, 1, 1], power = mean(ω.ld_E.power[h, :, :], dims=[2,3]))
    ld_H = (t = ω.ld_H.t[h, 1, 1], power = mean(ω.ld_H.power[h, :, :], dims=[2,3]))
    # Production
    pv = (t = ω.pv.t[h, 1, 1], power =  mean(ω.pv.power[h, :, :], dims=[2,3]), cost =  ω.pv.cost[y, s])
    # Electricity tariff
    grid = (cost_in = mean(ω.grid.cost_in[h, y, s], dims=[2,3]), cost_out =  mean(ω.grid.cost_out[h, y, s], dims=[2,3]))

    # Investment costs
    liion = (cost =  ω.liion.cost[y, s],)
    tes = (cost =  ω.tes.cost[y, s],)
    h2tank = (cost =  ω.h2tank.cost[y, s],)
    elyz = (cost =  ω.elyz.cost[y, s],)
    fc = (cost =  ω.fc.cost[y, s],)
    heater = (cost =  ω.heater.cost[y, s],)

    return Scenarios(ld_E, ld_H, pv, liion, tes, h2tank, elyz, fc, heater, grid)
end

function scenarios_reduction_clustering(ω::Scenarios, y::Int64, s::Int64, ncluster::Int64, algo::String)
    # Parmameters
    nh, ny, ns = size(ω.ld_E.power)

    # Clustering aggregated data
    results = clustering(ω.pv.power, ω.ld_E.power, ω.ld_H.power, ω.grid.cost_in, ω.grid.cost_out, ncluster = ncluster, algo = algo)
    # TODO rajouter proba
    if algo == "kmeans"
        # Demand
        ld_E = (t = repeat(ω.ld_E.t[:, 1, 1], 1, ncluster), power = results.centers[nh+1:2*nh,:])
        ld_H = (t = repeat(ω.ld_E.t[:, 1, 1], 1, ncluster), power = results.centers[2*nh+1:3*nh,:])
        # Production
        pv = (t = repeat(ω.ld_E.t[:, 1, 1], 1, ncluster), power = results.centers[1:nh,:], cost =  ω.pv.cost[y, s])
        # Electricity tariff
        grid = (cost_in = results.centers[3*nh+1:4*nh,:], cost_out = results.centers[4*nh+1:5*nh,:])
    elseif algo == "kmedoids"
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

    return Scenarios(ld_E, ld_H, pv, liion, tes, h2tank, elyz, fc, heater, grid)
end

# Scenario concatenation
function concatenate(ω1::Scenarios, ω2::Scenarios; dims::Int64 = 2)
    @assert 1 < dims < 4 "concatenation along y and s only !"
    # Demand
    ld_E = (t = cat(ω1.ld_E.t, ω2.ld_E.t, dims=dims), power = cat(ω1.ld_E.power, ω2.ld_E.power, dims=dims))
    ld_H = (t = cat(ω1.ld_H.t, ω2.ld_H.t, dims=dims), power = cat(ω1.ld_H.power, ω2.ld_H.power, dims=dims))
    # Production
    pv = (t = cat(ω1.pv.t, ω2.pv.t, dims=dims), power =  cat(ω1.pv.power, ω2.pv.power, dims=dims), cost =  cat(ω1.pv.cost, ω2.pv.cost, dims=dims-1))
    # Investment costs
    liion = (cost =  cat(ω1.liion.cost, ω2.liion.cost, dims=dims-1),)
    tes = (cost =  cat(ω1.tes.cost, ω2.tes.cost, dims=dims-1),)
    h2tank = (cost =  cat(ω1.h2tank.cost, ω2.h2tank.cost, dims=dims-1),)
    elyz = (cost =  cat(ω1.elyz.cost, ω2.elyz.cost, dims=dims-1),)
    fc = (cost =  cat(ω1.fc.cost, ω2.fc.cost, dims=dims-1),)
    heater = (cost =  cat(ω1.heater.cost, ω2.heater.cost, dims=dims-1),)
    # Electricity tariff
    grid = (cost_in =  cat(ω1.grid.cost_in, ω2.grid.cost_in, dims=dims), cost_out =  cat(ω1.grid.cost_out, ω2.grid.cost_out, dims=dims))

    return Scenarios(ld_E, ld_H, pv, liion, tes, h2tank, elyz, fc, heater, grid)
end
