#=
    Dummy designer
=#

mutable struct DummyDesigner <: AbstractDesigner
    u::NamedTuple
    DummyDesigner() = new()
end

### Offline
function initialize_designer!(des::DistributedEnergySystem, designer::DummyDesigner, ω::AbstractScenarios)
    # Preallocation
    preallocate!(designer, des.parameters.ny, des.parameters.ns)
    return designer
end

### Online
function compute_investment_decisions!(y::Int64, s::Int64, des::DistributedEnergySystem, designer::DummyDesigner)
    ϵ = 0.1

    if y != 1
        isa(des.liion, Liion) && des.liion.soh[end,y,s] < ϵ ? designer.u.liion[y,s] = designer.u.liion[1,s] : nothing
        isa(des.elyz, Electrolyzer) && des.elyz.soh[end,y,s] < ϵ ? designer.u.elyz[y,s] = designer.u.elyz[1,s]  : nothing
        isa(des.fc, FuelCell) && des.fc.soh[end,y,s] < ϵ ? designer.u.fc[y,s] = designer.u.fc[1,s]  : nothing
    end
end
