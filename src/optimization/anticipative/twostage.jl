#=
    Anticipative one-stage designer with anticipative controller
=#

mutable struct AnticipativeTwoStage <: AbstractDesigner
    options::MILPOptions
    u::NamedTuple
    model::JuMP.Model

    AnticipativeTwoStage(; options = MILPOptions()) = new(options)
end

### Offline
function offline_optimization!(des::DistributedEnergySystem, designer::AnticipativeTwoStage, ω::Scenarios)
    # Build anticipative controller
    controller = Anticipative()

    # Pre-allocate
    preallocate!(controller, des.parameters.nh, des.parameters.ny, des.parameters.ns)
    preallocate!(designer, des.parameters.ny, des.parameters.ns)

    # Scenario reduction from the optimization scenario pool
    println("Starting scenario reduction...")
    ω_reduced, probabilities = reduce(designer.options.reducer, ω)

    # Initialize model
    println("Building the model...")
    designer.model = build_model(des, MILP(options = designer.options), ω_reduced, probabilities)

    # Compute both investment and operation decisions
    println("Starting optimization...")
    optimize!(designer.model)

    # Assign controller values
    controller.u.liion .= value.(designer.model[:p_liion_ch] + designer.model[:p_liion_dch])
    controller.u.h2tank .= value.(designer.model[:p_h2tank_ch] + designer.model[:p_h2tank_dch])
    controller.u.tes .= value.(designer.model[:p_tes_ch] + designer.model[:p_tes_dch])
    controller.u.elyz .= value.(designer.model[:p_elyz_E])
    controller.u.fc .= value.(designer.model[:p_fc_E])
    controller.u.heater .= value.(designer.model[:p_heater_E])

    # Assign designer values
    designer.u.pv[1,:] .= value(designer.model[:r_pv])
    designer.u.liion[1,:] .= value(designer.model[:r_liion])
    designer.u.h2tank[1,:] .= value(designer.model[:r_h2tank])
    designer.u.elyz[1,:] .= value(designer.model[:r_elyz])
    designer.u.fc[1,:] .= value(designer.model[:r_fc])
    designer.u.tes[1,:] .= value(designer.model[:r_tes])

    return controller, designer
end

### Online
function compute_investment_decisions!(y::Int64, s::Int64, des::DistributedEnergySystem, designer::AnticipativeTwoStage)
    ϵ = 0.1

    isa(des.liion, Liion) && des.liion.soh[end,y,s] < ϵ ? designer.u.liion[y,s] = designer.u.liion[1,s] : nothing
    isa(des.elyz, Electrolyzer) && des.elyz.soh[end,y,s] < ϵ ? designer.u.elyz[y,s] = designer.u.elyz[1,s] : nothing
    isa(des.fc, FuelCell) && des.fc.soh[end,y,s] < ϵ ? designer.u.fc[y,s] = designer.u.fc[1,s] : nothing
end
