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
    isa(des.h2tank, H2Tank) ? controller.u.h2tank .= value.(designer.model[:p_h2tank_ch] + designer.model[:p_h2tank_dch]) : nothing
    isa(des.tes, ThermalSto) ? controller.u.tes .= value.(designer.model[:p_tes_ch] + designer.model[:p_tes_dch]) : nothing
    isa(des.elyz, Electrolyzer) ? controller.u.elyz .= value.(designer.model[:p_elyz_E]) : nothing
    isa(des.fc, FuelCell) ? controller.u.fc .= value.(designer.model[:p_fc_E]) : nothing
    isa(des.heater, Heater) ? controller.u.heater .= value.(designer.model[:p_heater_E]) : nothing

    # Preallocate and assigned values to the designer
    preallocate!(designer, des.parameters.ny, des.parameters.ns)
    isa(des.pv, Source) ? designer.u.pv[1,:] .= value(designer.model[:r_pv]) : nothing
    isa(des.liion, Liion) ? designer.u.liion[1,:] .= value(designer.model[:r_liion]) : nothing
    isa(des.h2tank, H2Tank) ? designer.u.h2tank[1,:] .= value(designer.model[:r_h2tank]) : nothing
    isa(des.elyz, Electrolyzer) ? designer.u.elyz[1,:] .= value(designer.model[:r_elyz]) : nothing
    isa(des.fc, FuelCell) ? designer.u.fc[1,:] .= value(designer.model[:r_fc]) : nothing
    isa(des.tes, ThermalSto) ? designer.u.tes[1,:] .= value(designer.model[:r_tes]) : nothing

    return controller, designer
end

### Online
function compute_investment_decisions!(y::Int64, s::Int64, des::DistributedEnergySystem, designer::AnticipativeEAC)
    ϵ = 0.1

    isa(des.liion, Liion) && des.liion.soh[end,y,s] < ϵ ? designer.u.liion[y,s] = des.liion.Erated[y,s] : nothing
    isa(des.elyz, Electrolyzer) && des.elyz.soh[end,y,s] < ϵ ? designer.u.elyz[y,s] = des.elyz.powerMax[y,s] : nothing
    isa(des.fc, FuelCell) && des.fc.soh[end,y,s] < ϵ ? designer.u.fc[y,s] = des.fc.powerMax[y,s] : nothing
end
