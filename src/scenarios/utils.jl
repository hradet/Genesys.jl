#=
    Utility functions
=#
# Concatenate scenario
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
# Repeat scenario
function Base.repeat(ω::Scenarios, nh::Int64, ny::Int64, ns::Int64)
    # Demand
    ld_E = (t = repeat(ω.ld_E.t, nh, ny, ns), power = repeat(ω.ld_E.power, nh, ny, ns))
    ld_H = (t = repeat(ω.ld_H.t, nh, ny, ns), power = repeat(ω.ld_H.power, nh, ny, ns))
    # Production
    pv = (t = repeat(ω.pv.t, nh, ny, ns), power =  repeat(ω.pv.power, nh, ny, ns), cost =  repeat(ω.pv.cost, ny, ns))
    # Investment costs
    liion = (cost =  repeat(ω.liion.cost, ny, ns),)
    tes = (cost =  repeat(ω.tes.cost, ny, ns),)
    h2tank = (cost =  repeat(ω.h2tank.cost, ny, ns),)
    elyz = (cost =  repeat(ω.elyz.cost, ny, ns),)
    fc = (cost =  repeat(ω.fc.cost, ny, ns),)
    heater = (cost =  repeat(ω.heater.cost, ny, ns),)
    # Electricity tariff
    grid = (cost_in =  repeat(ω.grid.cost_in, nh, ny, ns), cost_out =  repeat(ω.grid.cost_out, nh, ny, ns))

    return Scenarios(ld_E, ld_H, pv, liion, tes, h2tank, elyz, fc, heater, grid)
end
# Reshape scenario
function Base.reshape(ω::Scenarios, nh::Int64, ny::Int64, ns::Int64)
    ld_E = (t = reshape(ω.ld_E.t, nh, ny, ns), power = reshape(ω.ld_E.power, nh, ny, ns))
    ld_H = (t = reshape(ω.ld_H.t, nh, ny, ns), power = reshape(ω.ld_H.power, nh, ny, ns))
    # Production
    pv = (t = reshape(ω.pv.t, nh, ny, ns), power =  reshape(ω.pv.power, nh, ny, ns), cost = reshape(ω.pv.cost, ny, ns),)
    # Investment costs
    liion = (cost =  reshape(ω.liion.cost, ny, ns),)
    tes = (cost =  reshape(ω.tes.cost, ny, ns),)
    h2tank = (cost =  reshape(ω.h2tank.cost, ny, ns),)
    elyz = (cost =  reshape(ω.elyz.cost, ny, ns),)
    fc = (cost =  reshape(ω.fc.cost, ny, ns),)
    heater = (cost =  reshape(ω.heater.cost, ny, ns),)
    # Electricity tariff
    grid = (cost_in =  reshape(ω.grid.cost_in, nh, ny, ns), cost_out =  reshape(ω.grid.cost_out, nh, ny, ns))

    return Scenarios(ld_E, ld_H, pv, liion, tes, h2tank, elyz, fc, heater, grid)
end
# True if t is a weeken day
isweekend(t::Union{DateTime, Array{DateTime}}) = (Dates.dayname.(t) .== "Saturday") .| (Dates.dayname.(t) .== "Sunday")
# Chose the right markovchain as a function of t
chose(generator::MarkovGenerator, t::DateTime) = isweekend(t) ? generator.markovchain.wkd : generator.markovchain.wk
# Generate a yearly profile from typical days clustered data
generate(data_td::Array{Float64,2}, assignments::Array{Int64,1}) = reshape(data_td[:, assignments, :], size(data_td, 1) * size(assignments, 1))
