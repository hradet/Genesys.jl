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
