abstract type AbstractDemand end
abstract type AbstractGeneration end
abstract type AbstractStorage end
abstract type AbstractConverter end
abstract type AbstractGrid end
abstract type AbstractDesigner end
abstract type AbstractController end

mutable struct GlobalParameters
    ns::Int64 # number of [nh, ny] scenarios
    Δh::Int64 # operations time step in hours
    nh::Int64 # number of operation stages
    Δy::Int64 # investment time step in years
    ny::Int64 # number of investment stages
    τ::Float64 # discount rate
    renewable_share::Float64 # share of renewables [0,1]

    GlobalParameters(nh, ny, ns;
                Δh = 1,
                Δy = 1,
                τ = 0.045,
                renewable_share = 0.) = renewable_share == 1. ? new(ns, Δh, nh, Δy, ny, τ, 0.9999) : new(ns, Δh, nh, Δy, ny, τ, renewable_share)
end

mutable struct Microgrid
    parameters::GlobalParameters
    demands::Vector{AbstractDemand}
    generations::Vector{AbstractGeneration}
    storages::Vector{AbstractStorage}
    converters::Vector{AbstractConverter}
    grids::Vector{AbstractGrid}

    Microgrid(; parameters = GlobalParameters(8760, 20, 1)) = new(parameters)
end

# Add node to microgrid
# TODO : add asset by asset without initializing everything...
function add!(mg::Microgrid, assets...)
    # Build and preallocate
    mg.demands = [preallocate!(a, mg.parameters.nh, mg.parameters.ny, mg.parameters.ns) for a in assets if a isa AbstractDemand]
    mg.generations = [preallocate!(a, mg.parameters.nh, mg.parameters.ny, mg.parameters.ns) for a in assets if a isa AbstractGeneration]
    mg.storages = [preallocate!(a, mg.parameters.nh, mg.parameters.ny, mg.parameters.ns) for a in assets if a isa AbstractStorage]
    mg.converters = [preallocate!(a, mg.parameters.nh, mg.parameters.ny, mg.parameters.ns) for a in assets if a isa AbstractConverter]
    mg.grids = [preallocate!(a, mg.parameters.nh, mg.parameters.ny, mg.parameters.ns) for a in assets if a isa AbstractGrid]
end

# Preallocate abstract designer
function preallocate!(mg::Microgrid, designer::AbstractDesigner)
    designer.decisions = (generations = [zeros(mg.parameters.ny, mg.parameters.ns) for a in mg.generations],
                          storages = [zeros(mg.parameters.ny, mg.parameters.ns) for a in mg.storages],
                          converters = [zeros(mg.parameters.ny, mg.parameters.ns) for a in mg.converters])
end
# Preallocate abstract controller
function preallocate!(mg::Microgrid, controller::AbstractController)
    controller.decisions = (converters = [zeros(mg.parameters.nh, mg.parameters.ny, mg.parameters.ns) for a in mg.converters],
                            storages = [zeros(mg.parameters.nh, mg.parameters.ny, mg.parameters.ns) for a in mg.storages])
end

# Copy function
function Base.copy(mg::Microgrid, nh::Int64, ny::Int64, ns::Int64)
    # TODO : copy not working without default parameters !
    # Initialize mg
    microgrid = Microgrid(parameters = GlobalParameters(nh, ny, ns, renewable_share = mg.parameters.renewable_share))
    demands, generations, storages, converters, grids = [], [], [], [], []
    # Get the assets
    # Demands
    for a in mg.demands
        if a.carrier isa Electricity
            push!(demands, Demand(carrier = Electricity()))
        elseif a.carrier isa Heat
            push!(demands, Demand(carrier = Heat()))
        elseif a.carrier isa Hydrogen
            push!(demands, Demand(carrier = Hydrogen()))
        end
    end
    # Generations
    for a in mg.generations
        if a isa Solar
            push!(generations, Solar())
        end
    end
    # Storages
    for a in mg.storages
        if a isa Liion
            push!(storages, Liion())
        elseif a isa ThermalStorage
            push!(storages, ThermalStorage())
        elseif a isa H2Tank
            push!(storages, H2Tank())
        end
    end
    # Converters
    for a in mg.converters
        if a isa Heater
            push!(converters, Heater())
        elseif a isa Electrolyzer
            push!(converters, Electrolyzer())
        elseif a isa FuelCell
            push!(converters, FuelCell())
        end
    end
    # Grids
    for a in mg.grids
        if a.carrier isa Electricity
            push!(grids, Grid(carrier = Electricity()))
        elseif a.carrier isa Heat
            push!(grids, Grid(carrier = Heat()))
        elseif a.carrier isa Hydrogen
            push!(grids, Grid(carrier = Hydrogen()))
        end
    end
    # Build the microgrid
    add!(microgrid, demands..., generations..., storages..., converters..., grids...)

    return microgrid
end

# Find if the datatype is in a mg field
function isin(field::Vector, type::DataType)
    # Return true if the datatype is in the field and its index
    bool, idx = false, NaN
    for (k, a) in enumerate(field)
        if a isa type
            bool, idx = true, k
        end
    end
    return bool, idx
end
