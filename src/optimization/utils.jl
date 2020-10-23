abstract type AbstractOneStageDesigner <: AbstractDesigner end
abstract type AbstractOneStageStochasticDesigner <: AbstractDesigner end
abstract type AbstractMultiStageDesigner <: AbstractDesigner end
abstract type AbstractMultiStageStochasticDesigner <: AbstractDesigner end

function set_bounds(des::DistributedEnergySystem)
    lb, ub = Float64[], Float64[]
    if isa(des.pv, Source)
        push!(lb, 0.) ; push!(ub, 1000.)
    end
    if isa(des.liion, Liion)
        push!(lb, 0.) ; push!(ub, 1000.)
    end
    if isa(des.h2tank, H2Tank)
        push!(lb, 0.) ; push!(ub, 50000.)
    end
    if isa(des.elyz, Electrolyzer)
        push!(lb, 0.) ; push!(ub, 50.)
    end
    if isa(des.fc, FuelCell)
        push!(lb, 0.) ; push!(ub, 50.)
    end
    if isa(des.tes, ThermalSto)
        push!(lb, 0.) ; push!(ub, 1000.)
    end

    return lb, ub
end
