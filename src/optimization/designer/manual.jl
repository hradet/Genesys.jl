#=
    Manual designer
=#

mutable struct Manual <: AbstractDesigner
    generations::Vector{Float64}
    storages::Vector{Float64}
    converters::Vector{Float64}
    decisions::NamedTuple

    Manual(; generations = [0.], storages = [0.], converters = [0.]) = new(generations, storages, converters)
end

### Offline
function initialize_designer!(mg::Microgrid, designer::Manual, ω::AbstractScenarios)
    # Preallocation
    preallocate!(mg, designer)

    # Fix initial values
    for (k, a) in enumerate(designer.decisions.generations)
        a[1,:] .= designer.generations[k]
    end
    for (k, a) in enumerate(designer.decisions.storages)
        a[1,:] .= designer.storages[k]
    end
    for (k, a) in enumerate(designer.decisions.converters)
        a[1,:] .= designer.converters[k]
    end

    return designer
end

### Online
function compute_investment_decisions!(y::Int64, s::Int64, mg::Microgrid, designer::Manual)
    # TODO
    return nothing
end
