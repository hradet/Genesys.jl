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
