#=
    Scenario reduction functions
=#
abstract type AbstractScenarios end

mutable struct Scenarios <: AbstractScenarios
    demands::Vector{Any}
    generations::Vector{Any}
    storages::Vector{Any}
    converters::Vector{Any}
    grids::Vector{Any}
end

# Constructor
function Scenarios(d::Dict{}, mg::Microgrid)
    # Utils to simplify the writting
    h, y, s = 1:mg.parameters.nh, 1:mg.parameters.ny, 1:mg.parameters.ns
    # Initialize
    demands, generations, storages, converters, grids = [], [], [], [], []
    # Demands
    for (k, a) in enumerate(mg.demands)
        if a.carrier isa Electricity
            push!(demands, (t = d["ld_E"].t[h, y, s], power = d["ld_E"].power[h, y, s]))
        elseif a.carrier isa Heat
            push!(demands, (t = d["ld_H"].t[h, y, s], power = d["ld_H"].power[h, y, s]))
        end
    end
    # Generation
    for (k, a) in enumerate(mg.generations)
        if a isa Solar
            push!(generations, (t = d["pv"].t[h, y, s], power = d["pv"].power[h, y, s], cost = d["pv"].cost[y, s]))
        end
    end
    # Storages
    for (k, a) in enumerate(mg.storages)
        if a isa Liion
            push!(storages, (cost = d["liion"].cost[y, s],))
        elseif a isa ThermalStorage
            push!(storages, (cost = d["tes"].cost[y, s],))
        elseif a isa H2Tank
            push!(storages, (cost = d["h2tank"].cost[y, s],))
        end
    end
    # Converters
    for (k, a) in enumerate(mg.converters)
        if a isa Electrolyzer
            push!(converters, (cost = d["elyz"].cost[y, s],))
        elseif a isa FuelCell
            push!(converters, (cost = d["fc"].cost[y, s],))
        elseif a isa Heater
            push!(converters, (cost = d["heater"].cost[y, s],))
        end
    end
    # Grids
    for (k, a) in enumerate(mg.grids)
        if a.carrier isa Electricity
            push!(grids, (cost_in = d["grid"].cost_in[h, y, s], cost_out = d["grid"].cost_out[h, y, s]))
        end
    end

    return Scenarios(demands, generations, storages, converters, grids)
end
