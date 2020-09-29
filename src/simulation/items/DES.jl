# Distributed energy system
mutable struct GlobalParameters
    ns::Int64 # number of [nh, ny] scenarios
    Δh::Int64 # operations time step in hours
    nh::Int64 # number of operation stages
    Δy::Int64 # investment time step in years
    ny::Int64 # number of investment stages
    τ::Float64 # discount rate
    τ_share::Float64 # share of renewables [0,1]

    GlobalParameters(;ns = 1,
                Δh = 1,
                nh = 8760,
                Δy = 1,
                ny = 20,
                τ = 0.045,
                τ_share = 0.) =
                new(ns, Δh, nh, Δy, ny, τ, τ_share)
end

mutable struct DistributedEnergySystem
    parameters::GlobalParameters
    # Loads
    ld_E::Union{Nothing, Load}
    ld_H::Union{Nothing, Load}
    # Source
    pv::Union{Nothing, Source}
    # Storage
    liion::Union{Nothing, Liion}
    tes::Union{Nothing, ThermalSto}
    h2tank::Union{Nothing, H2Tank}
    # Converters
    heater::Union{Nothing, Heater}
    elyz::Union{Nothing, Electrolyzer}
    fc::Union{Nothing, FuelCell}
    # Grid
    grid::Union{Nothing, Grid}
end

function DistributedEnergySystem(; ld_E = nothing,
                                   ld_H = nothing,
                                   pv = nothing,
                                   liion = nothing,
                                   tes = nothing,
                                   h2tank = nothing,
                                   heater = nothing,
                                   elyz = nothing,
                                   fc = nothing,
                                   grid = nothing,
                                   parameters = GlobalParameters())
    # Parameters
    nh = parameters.nh
    ny = parameters.ny
    ns = parameters.ns

    # Loads
    isa(ld_E, Load) ? preallocate!(ld_E, nh, ny, ns) : nothing
    isa(ld_H, Load) ? preallocate!(ld_H, nh, ny, ns) : nothing
    # Source
    isa(pv, Source) ? preallocate!(pv, nh, ny, ns) : nothing
    # Storage
    isa(liion, Liion) ? preallocate!(liion, nh, ny, ns) : nothing
    isa(tes, ThermalSto) ? preallocate!(tes, nh, ny, ns) : nothing
    isa(h2tank, H2Tank) ? preallocate!(h2tank, nh, ny, ns) : nothing
    # Converters
    isa(heater, Heater) ? preallocate!(heater, nh, ny, ns) : nothing
    isa(elyz, Electrolyzer) ? preallocate!(elyz, nh, ny, ns) : nothing
    isa(fc, FuelCell) ? preallocate!(fc, nh, ny, ns) : nothing
    # Grid
    isa(grid, Grid) ? preallocate!(grid, nh, ny, ns) : nothing

    return DistributedEnergySystem(parameters, ld_E, ld_H, pv, liion, tes, h2tank, heater, elyz, fc, grid)
end
