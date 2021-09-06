#=
    Utility functions
=#
# Concatenate scenario
function Base.cat(ω1::Scenarios, ω2::Scenarios; dims::Int64 = 2)
    @assert 1 < dims < 4 "concatenation along dims 2 and 3 only !"
    # Initialize
    demands, generations, storages, converters, grids = [], [], [], [], []
    # Demands
    for k in length(ω1.demands)
        push!(demands, (t = cat(ω1.demands[k].t, ω2.demands[k].t, dims=dims), power = cat(ω1.demands[k].power,  ω2.demands[k].power, dims=dims)))
    end
    # Generations
    for k in length(ω1.generations)
        push!(generations, (t = cat(ω1.generations[k].t, ω2.generations[k].t, dims=dims), power =  cat(ω1.generations[k].power, ω2.generations[k].power, dims=dims), cost =  cat(ω1.generations[k].cost, ω2.generations[k].cost, dims=dims-1)))
    end
    # Storages
    for k in length(ω1.storages)
        push!(storages, (cost =  cat(ω1.storages[k].cost, ω2.storages[k].cost, dims=dims-1),))
    end
    # Converters
    for k in length(ω1.converters)
        push!(converters, (cost =  cat(ω1.converters[k].cost, ω2.converters[k].cost, dims=dims-1),))
    end
    # Grids
    for k in length(ω1.grids)
        push!(grids, (cost_in = cat(ω1.grids[k].cost_in, ω2.grids[k].cost_in, dims=dims), cost_out =  cat(ω1.grids[k].cost_out, ω2.grids[k].cost_out, dims=dims)))
    end

    return Scenarios(demands, generations, storages, converters, grids)
end
# Repeat scenario
function Base.repeat(ω::Scenarios, nh::Int64, ny::Int64, ns::Int64)
    # Initialize
    demands, generations, storages, converters, grids = [], [], [], [], []
    # Demands
    for k in length(ω.demands)
        push!(demands, (t = repeat(ω.demands[k].t, nh, ny, ns), power = repeat(ω.demands[k].power, nh, ny, ns)))
    end
    # Generations
    for k in length(ω.generations)
        push!(generations, (t = repeat(ω.generations[k].t, nh, ny, ns), power =  repeat(ω.generations[k].power, nh, ny, ns), cost =  repeat(ω.generations[k].cost, ny, ns)))
    end
    # Storages
    for k in length(ω.storages)
        push!(storages, (cost =  repeat(ω.storages[k].cost, ny, ns),))
    end
    # Converters
    for k in length(ω.converters)
        push!(converters, (cost =  repeat(ω.converters[k].cost, ny, ns),))
    end
    # Grids
    for k in length(ω.grids)
        push!(grids, (cost_in =  repeat(ω.grids[k].cost_in, nh, ny, ns), cost_out =  repeat(ω.grids[k].cost_out, nh, ny, ns)))
    end

    return Scenarios(demands, generations, storages, converters, grids)
end
# Reshape scenario
function Base.reshape(ω::Scenarios, nh::Int64, ny::Int64, ns::Int64)
    # Initialize
    demands, generations, storages, converters, grids = [], [], [], [], []
    # Demands
    for k in length(ω.demands)
        push!(demands, (t = reshape(ω.demands[k].t, nh, ny, ns), power = reshape(ω.demands[k].power, nh, ny, ns)))
    end
    # Generations
    for k in length(ω.generations)
        push!(generations, (t = reshape(ω.generations[k].t, nh, ny, ns), power =  reshape(ω.generations[k].power, nh, ny, ns), cost =  reshape(ω.generations[k].cost, ny, ns)))
    end
    # Storages
    for k in length(ω.storages)
        push!(storages, (cost =  reshape(ω.storages[k].cost, ny, ns),))
    end
    # Converters
    for k in length(ω.converters)
        push!(converters, (cost =  reshape(ω.converters[k].cost, ny, ns),))
    end
    # Grids
    for k in length(ω.grids)
        push!(grids, (cost_in =  reshape(ω.grids[k].cost_in, nh, ny, ns), cost_out =  reshape(ω.grids[k].cost_out, nh, ny, ns)))
    end

    return Scenarios(demands, generations, storages, converters, grids)
end
# True if t is a weeken day
isweekend(t::Union{DateTime, Array{DateTime}}) = (Dates.dayname.(t) .== "Saturday") .| (Dates.dayname.(t) .== "Sunday")
# Chose the right markovchain as a function of t
chose(generator::MarkovGenerator, t::DateTime) = isweekend(t) ? generator.markovchain.wkd : generator.markovchain.wk
# Generate a yearly profile from typical days clustered data
generate(data_td::Array{Float64,2}, assignments::Array{Int64,1}) = reshape(data_td[:, assignments, :], size(data_td, 1) * size(assignments, 1))
# Normalization
min_max_normalization(x::AbstractArray) = maximum(x) == minimum(x) ? x ./ maximum(x) : (x .- minimum(x)) ./ (maximum(x) .- minimum(x))
