# Distributed energy system
mutable struct GlobalParameters
    ns::Int64 # number of [nh, ny] scenarios
    Δh::Int64 # Operations time step in hours
    H::Int64 # Operation horizon in hours
    Δy::Int64 # Investment time step in years
    Y::Int64 # Investment horizon in years
    τ::Float64 # Discount rate
    τ_share::Float64 # Share of renewables [0,1]

    GlobalParameters(;ns = 1,
                Δh = 1,
                H = 8760,
                Δy = 1,
                Y = 20,
                τ = 0.045,
                τ_share = 0.) =
                new(ns, Δh, H, Δy, Y, τ, τ_share)
end

mutable struct DES
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
    # Controller
    controller::Union{Nothing, AbstractController}
    designer::Union{Nothing, AbstractDesigner}
end

function DES(; ld_E = nothing,
            ld_H = nothing,
            pv = nothing,
            liion = nothing,
            tes = nothing,
            h2tank = nothing,
            heater = nothing,
            elyz = nothing,
            fc = nothing,
            grid = nothing,
            controller = DummyController(),
            designer = DummyDesigner(),
            parameters = GlobalParameters())
    # Parameters
    nh = length(parameters.Δh : parameters.Δh : parameters.H)
    ny = length(parameters.Δy : parameters.Δy : parameters.Y)
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
    # Controller
    preallocate!(controller, nh, ny, ns)
    # Designer
    preallocate!(designer, nh, ny, ns)

    return DES(parameters, ld_E, ld_H, pv, liion, tes, h2tank, heater, elyz, fc, grid, controller, designer)
end
