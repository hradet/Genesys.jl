abstract type AbstractOneStageDesigner <: AbstractDesigner end
abstract type AbstractOneStageStochasticDesigner <: AbstractDesigner end
abstract type AbstractMultiStageDesigner <: AbstractDesigner end
abstract type AbstractMultiStageStochasticDesigner <: AbstractDesigner end

function set_bounds(des::DistributedEnergySystem)

    lb, ub = zeros(6), zeros(6)
    isa(des.pv, Source) ? ub[1] = 1000. : nothing
    isa(des.liion, Liion) ? ub[2] = 1000. : nothing
    isa(des.h2tank, H2Tank) ? ub[3] = 50000. : nothing
    isa(des.elyz, Electrolyzer) ? ub[4] = 50. : nothing
    isa(des.fc, FuelCell) ? ub[5] = 50. : nothing
    isa(des.tes, ThermalSto) ? ub[6] = 1000. : nothing

    return lb, ub
end

function copy(des::DistributedEnergySystem, nh::Int64, ny::Int64, ns::Int64)
    des_copy = DistributedEnergySystem(ld_E = isa(des.ld_E, Load) ? Load() : nothing,
                                  ld_H = isa(des.ld_H, Load) ? Load() : nothing,
                                  pv = isa(des.pv, Source) ? Source() : nothing,
                                  liion = isa(des.liion, Liion) ? Liion() : nothing,
                                  tes = isa(des.tes, ThermalSto) ? ThermalSto() : nothing,
                                  h2tank = isa(des.h2tank, H2Tank) ? H2Tank() : nothing,
                                  elyz = isa(des.elyz, Electrolyzer) ? Electrolyzer() : nothing,
                                  fc = isa(des.fc, FuelCell) ? FuelCell() : nothing,
                                  heater = isa(des.heater, Heater) ? Heater() : nothing,
                                  grid = isa(des.grid, Grid) ? Grid() : nothing,
                                  parameters = Genesys.GlobalParameters(nh, ny, ns, τ_share = des.parameters.τ_share))

    return des_copy
end
