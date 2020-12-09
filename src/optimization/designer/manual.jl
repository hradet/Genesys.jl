#=
    Manual designer
=#

mutable struct Manual <: AbstractDesigner
    pv::Float64
    liion::Float64
    tes::Float64
    h2tank::Float64
    elyz::Float64
    fc::Float64
    u::NamedTuple

    Manual(; pv = 0.,
             liion = 0.,
             tes = 0.,
             h2tank = 0.,
             elyz = 0.,
             fc = 0.) =
             new(pv, liion, tes, h2tank, elyz, fc)
end

### Offline
function initialize_designer!(des::DistributedEnergySystem, designer::Manual, ω::AbstractScenarios)
    # Preallocation
    preallocate!(designer, des.parameters.ny, des.parameters.ns)

    # Fix initial values
    designer.u.pv[1,:] .= designer.pv
    designer.u.liion[1,:] .= designer.liion
    designer.u.tes[1,:] .= designer.tes
    designer.u.h2tank[1,:] .= designer.h2tank
    designer.u.elyz[1,:] .= designer.elyz
    designer.u.fc[1,:] .= designer.fc

    return designer
end

### Online
function compute_investment_decisions!(y::Int64, s::Int64, des::DistributedEnergySystem, designer::Manual)
    ϵ = 0.1

    if y != 1
        isa(des.liion, Liion) && des.liion.soh[end,y,s] < ϵ ? designer.u.liion[y,s] = designer.u.liion[1,s] : nothing
        isa(des.elyz, Electrolyzer) && des.elyz.soh[end,y,s] < ϵ ? designer.u.elyz[y,s] = designer.u.elyz[1,s]  : nothing
        isa(des.fc, FuelCell) && des.fc.soh[end,y,s] < ϵ ? designer.u.fc[y,s] = designer.u.fc[1,s]  : nothing
    end
end
