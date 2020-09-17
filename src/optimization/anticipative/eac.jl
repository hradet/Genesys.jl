#=
    Anticipative one-stage designer with anticipative controller
=#

mutable struct AnticipativeEAC <: AbstractOneStageStochasticDesigner
    options::EACStochOptions
    u::NamedTuple
    model::JuMP.Model
    AnticipativeEAC(; options = EACStochOptions()) = new(options)
end

### Offline
function offline_optimization!(des::DistributedEnergySystem, designer::AnticipativeEAC, ω::AbstractScenarios)
    # Scenario reduction from the optimization scenario pool
    ω_anticipative = scenarios_reduction(EACStoch(options = designer.options), ω)

    # Build model
    designer.model = build_model(des, EACStoch(options = designer.options), ω_anticipative)

    # Compute both investment and operation decisions
    optimize!(designer.model)

    # Preallocate and assigned values to the controller
    controller = Anticipative()
    preallocate!(controller, des.parameters.nh, des.parameters.ny, des.parameters.ns)
    isa(des.liion, Liion) ? controller.u.liion .= value.(designer.model[:p_liion_ch] + designer.model[:p_liion_dch]) : nothing

    # Preallocate and assigned values to the designer
    preallocate!(designer, des.parameters.ny, des.parameters.ns)
    isa(des.liion, Liion) ? designer.u.liion[1] = value.(designer.model[:r_liion]) : nothing
    isa(des.pv, Source) ? designer.u.pv[1] = value.(designer.model[:r_pv]) : nothing

    return controller, designer
end

### Online
function compute_investment_decisions!(y::Int64, s::Int64, des::DistributedEnergySystem, designer::AnticipativeEAC)
    ϵ = 0.1

    # Liion
    if isa(des.liion, Liion) && des.liion.soh[end,y,s] < ϵ
        designer.u.liion[y,s] = des.liion.Erated[y,s]
    end

    # Electrolyzer
    if isa(des.elyz, Electrolyzer) && des.elyz.soh[end,y,s] < ϵ
        designer.u.elyz[y,s] = des.elyz.powerMax[y,s]
    end

    # FuelCell
    if isa(des.fc, FuelCell) && des.fc.soh[end,y,s] < ϵ
        designer.u.fc[y,s] = des.fc.powerMax[y,s]
    end
end
