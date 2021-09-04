
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
                renewable_share = 0.) =
                new(ns, Δh, nh, Δy, ny, τ, renewable_share)
end

mutable struct Microgrid
    parameters::GlobalParameters
    demands::Vector{Any}
    generations::Vector{Any}
    storages::Vector{Any}
    converters::Vector{Any}
    grids::Vector{Any}

    Microgrid(; parameters = GlobalParameters(8760, 20, 1)) = new(parameters)
end

# Add node to microgrid
function add!(mg::Microgrid, assets...)
    # Build and preallocate
    mg.demands = [preallocate!(a, mg.parameters.nh, mg.parameters.ny, mg.parameters.ns) for a in assets if a isa AbstractDemand]
    mg.generations = [preallocate!(a, mg.parameters.nh, mg.parameters.ny, mg.parameters.ns) for a in assets if a isa AbstractGeneration]
    mg.storages = [preallocate!(a, mg.parameters.nh, mg.parameters.ny, mg.parameters.ns) for a in assets if a isa AbstractStorage]
    mg.converters = [preallocate!(a, mg.parameters.nh, mg.parameters.ny, mg.parameters.ns) for a in assets if a isa AbstractConverter]
    mg.grids = [preallocate!(a, mg.parameters.nh, mg.parameters.ny, mg.parameters.ns) for a in assets if a isa AbstractGrid]
end


#
# function Base.copy(des::DistributedEnergySystem, nh::Int64, ny::Int64, ns::Int64)
#     des_copy = DistributedEnergySystem(ld_E = isa(des.ld_E, Load) ? Load() : nothing,
#                                   ld_H = isa(des.ld_H, Load) ? Load() : nothing,
#                                   pv = isa(des.pv, Source) ? Source() : nothing,
#                                   liion = isa(des.liion, Liion) ? Liion() : nothing,
#                                   tes = isa(des.tes, ThermalSto) ? ThermalSto() : nothing,
#                                   h2tank = isa(des.h2tank, H2Tank) ? H2Tank() : nothing,
#                                   elyz = isa(des.elyz, Electrolyzer) ? Electrolyzer() : nothing,
#                                   fc = isa(des.fc, FuelCell) ? FuelCell() : nothing,
#                                   heater = isa(des.heater, Heater) ? Heater() : nothing,
#                                   grid = isa(des.grid, Grid) ? Grid() : nothing,
#                                   parameters = Genesys.GlobalParameters(nh, ny, ns, renewable_share = des.parameters.renewable_share))
#
#     return des_copy
# end
