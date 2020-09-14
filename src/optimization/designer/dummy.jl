#=
    Dummy designer
=#

mutable struct DummyDesigner <: AbstractDesigner
    u::NamedTuple
    DummyDesigner() = new()
end

### Offline
function initialize_designer!(des::DES, ω::Scenarios)
    return nothing
end

### Online
function compute_investment_decisions!(y::Int64, s::Int64, des::DES)
    ϵ = 0.1

    # Liion
    if isa(des.liion, Liion) && des.liion.soh[end,y,s] < ϵ
        des.designer.u.liion[y,s] = des.liion.Erated[y,s]
    end

    # Electrolyzer
    if isa(des.elyz, Electrolyzer) && des.elyz.soh[end,y,s] < ϵ
        des.designer.u.elyz[y,s] = des.elyz.powerMax[y,s]
    end

    # FuelCell
    if isa(des.fc, FuelCell) && des.fc.soh[end,y,s] < ϵ
        des.designer.u.fc[y,s] = des.fc.powerMax[y,s]
    end
end
