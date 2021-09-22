#=
    Scenario reduction functions
=#
abstract type AbstractScenarios end

mutable struct Scenarios{T, O, I} <: AbstractScenarios
    demands::Vector{NamedTuple{(:t, :power),Tuple{T,O}}}
    generations::Vector{NamedTuple{(:t, :power, :cost), Tuple{T, O, I}}}
    storages::Vector{NamedTuple{(:cost,), Tuple{I}}}
    converters::Vector{NamedTuple{(:cost,), Tuple{I}}}
    grids::Vector{NamedTuple{(:cost_in, :cost_out), Tuple{O, O}}}
end

# Constructor
function Scenarios(d::Dict{}, mg::Microgrid)
    # Utils to simplify the writting
    h, y, s = 1:mg.parameters.nh, 1:mg.parameters.ny, 1:mg.parameters.ns
    T, O, I = Array{DateTime,3}, Array{Float64, 3}, Array{Float64, 2}
    # Initialize
    demands = Vector{NamedTuple{(:t, :power),Tuple{T,O}}}(undef, length(mg.demands))
    generations = Vector{NamedTuple{(:t, :power, :cost), Tuple{T, O, I}}}(undef, length(mg.generations))
    storages = Vector{NamedTuple{(:cost,), Tuple{I}}}(undef, length(mg.storages))
    converters = Vector{NamedTuple{(:cost,), Tuple{I}}}(undef, length(mg.converters))
    grids = Vector{NamedTuple{(:cost_in, :cost_out), Tuple{O, O}}}(undef, length(mg.grids))
    # Demands
    for (k, a) in enumerate(mg.demands)
        if a.carrier isa Electricity
            demands[k] = (t = d["ld_E"].t[h, y, s], power = d["ld_E"].power[h, y, s])
        elseif a.carrier isa Heat
            demands[k] = (t = d["ld_H"].t[h, y, s], power = d["ld_H"].power[h, y, s])
        end
    end
    # Generation
    for (k, a) in enumerate(mg.generations)
        if a isa Solar
            generations[k] = (t = d["pv"].t[h, y, s], power = d["pv"].power[h, y, s], cost = d["pv"].cost[y, s])
        end
    end
    # Storages
    for (k, a) in enumerate(mg.storages)
        if a isa Liion
            storages[k] = (cost = d["liion"].cost[y, s],)
        elseif a isa ThermalStorage
            storages[k] = (cost = d["tes"].cost[y, s],)
        elseif a isa H2Tank
            storages[k] = (cost = d["h2tank"].cost[y, s],)
        end
    end
    # Converters
    for (k, a) in enumerate(mg.converters)
        if a isa Electrolyzer
            converters[k] = (cost = d["elyz"].cost[y, s],)
        elseif a isa FuelCell
            converters[k] = (cost = d["fc"].cost[y, s],)
        elseif a isa Heater
            converters[k] = (cost = d["heater"].cost[y, s],)
        end
    end
    # Grids
    for (k, a) in enumerate(mg.grids)
        if a.carrier isa Electricity
            grids[k] = (cost_in = d["grid"].cost_in[h, y, s], cost_out = d["grid"].cost_out[h, y, s])
        end
    end

    return Scenarios(demands, generations, storages, converters, grids)
end
