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

# Scenario concatenation
function Base.cat(ω1::Scenarios, ω2::Scenarios; dims::Int64 = 2)
    @assert 1 < dims < 4 "concatenation along dims 2 and 3 only !"
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

function Base.reshape(ω::Scenarios{Array{DateTime,2}, Array{Float64,2}, Float64}, nh::Int64, ny::Int64, ns::Int64)
    ld_E = (t = reshape(ω.ld_E.t, nh, ny, ns), power = reshape(ω.ld_E.power, nh, ny, ns))
    ld_H = (t = reshape(ω.ld_H.t, nh, ny, ns), power = reshape(ω.ld_H.power, nh, ny, ns))
    # Production
    pv = (t = reshape(ω.pv.t, nh, ny, ns), power =  reshape(ω.pv.power, nh, ny, ns), cost = fill(ω.pv.cost, ny, ns),)
    # Investment costs
    liion = (cost =  fill(ω.liion.cost, ny, ns),)
    tes = (cost =  fill(ω.tes.cost, ny, ns),)
    h2tank = (cost =  fill(ω.h2tank.cost, ny, ns),)
    elyz = (cost =  fill(ω.elyz.cost, ny, ns),)
    fc = (cost =  fill(ω.fc.cost, ny, ns),)
    heater = (cost =  fill(ω.heater.cost, ny, ns),)
    # Electricity tariff
    grid = (cost_in =  reshape(ω.grid.cost_in, nh, ny, ns), cost_out =  reshape(ω.grid.cost_out, nh, ny, ns))

    return Scenarios(ld_E, ld_H, pv, liion, tes, h2tank, elyz, fc, heater, grid)
end
