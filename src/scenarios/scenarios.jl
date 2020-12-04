#=
    Scenario reduction functions
=#
abstract type AbstractScenarios end

mutable struct Scenarios{T, O, I} <: AbstractScenarios
    # Demand
    ld_E::NamedTuple{(:t, :power), Tuple{T, O}}
    ld_H::NamedTuple{(:t, :power), Tuple{T, O}}
    # Production
    pv::NamedTuple{(:t, :power, :cost), Tuple{T, O, I}}
    # Investment costs
    liion::NamedTuple{(:cost,), Tuple{I}}
    tes::NamedTuple{(:cost,), Tuple{I}}
    h2tank::NamedTuple{(:cost,), Tuple{I}}
    elyz::NamedTuple{(:cost,), Tuple{I}}
    fc::NamedTuple{(:cost,), Tuple{I}}
    heater::NamedTuple{(:cost,), Tuple{I}}
    # Electricity tariff
    grid::NamedTuple{(:cost_in, :cost_out), Tuple{O, O}}
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
