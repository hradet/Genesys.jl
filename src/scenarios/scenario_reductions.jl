#=
    Scenario reduction functions
=#
mutable struct Scenarios <: AbstractScenarios
     timestamp::Array{DateTime,1}
     values::NamedTuple
     proba::Array{Float64,1}
end

# Constructor
function Scenarios(outputGUI::Dict{}, nh::Int64, ny::Int64, ns::Int64)
    # Values
    values = (
    # Demand
    ld_E = outputGUI["ld_E"][1:nh, 1:ny, 1:ns],
    ld_H =outputGUI["ld_H"][1:nh, 1:ny, 1:ns],
    # Production
    pv_E = outputGUI["pv_E"][1:nh, 1:ny, 1:ns],
    # Investment costs
    C_pv = outputGUI["C_pv"][1:ny, 1:ns],
    C_liion = outputGUI["C_liion"][1:ny, 1:ns],
    C_tes = outputGUI["C_tes"][1:ny, 1:ns],
    C_tank = outputGUI["C_tank"][1:ny, 1:ns],
    C_elyz = outputGUI["C_elyz"][1:ny, 1:ns],
    C_fc = outputGUI["C_fc"][1:ny, 1:ns],
    C_heater = outputGUI["C_heater"][1:ny, 1:ns],
    # Electricity tariff
    C_grid_in = outputGUI["C_grid_in"][1:nh, 1:ny, 1:ns],
    C_grid_out = outputGUI["C_grid_out"][1:nh, 1:ny, 1:ns],
    )
    # Probabilities - uniform
    proba = 1 / ns * ones(ns)

    return Scenarios(outputGUI["timestamp"],values, proba)
end

# Scenario Reduction functions
# One stage designer reduction
function scenarios_reduction(designer::AbstractOneStageDesigner, ω::Scenarios)

    if designer.options.scenario_reduction == "manual"
        # Year and scenario indexes are manually chosen
        y =  designer.options.y
        s =  designer.options.s
    else
        # Year and scenario indexes are randomly chosen to 1
        y, s = 1, 1
    end

    # Reduced values
    values = (
    ld_E = ω.values.ld_E[:,y,s],
    ld_H = ω.values.ld_H[:,y,s],
    # Production
    pv_E = ω.values.pv_E[:,y,s],
    # Investment costs
    C_pv = ω.values.C_pv[y,s],
    C_liion = ω.values.C_liion[y,s],
    C_tes = ω.values.C_tes[y,s],
    C_tank = ω.values.C_tank[y,s],
    C_elyz = ω.values.C_elyz[y,s],
    C_fc = ω.values.C_fc[y,s],
    C_heater = ω.values.C_heater[y,s],
    # Electricity tariff
    C_grid_in = ω.values.C_grid_in[:,y,s],
    C_grid_out = ω.values.C_grid_out[:,y,s],
    )

    return Scenarios(ω.timestamp,values, [ω.proba[s]])
end
# One stage stochastic designer reduction
function scenarios_reduction(designer::AbstractOneStageStochasticDesigner, ω::Scenarios)

    if designer.options.scenario_reduction == "manual"
        # Year and scenario indexes are manually chosen
        y = designer.options.range_y
        s = designer.options.s
    else
        #TODO !! so far, year and scenario indexes are randomly chosen
        y, s = size(ω.values.ld_E,2), 1
    end

    # Reduced values
    values = (
    ld_E = ω.values.ld_E[:,y,s],
    ld_H = ω.values.ld_H[:,y,s],
    # Production
    pv_E = ω.values.pv_E[:,y,s],
    # Investment costs
    C_pv = ω.values.C_pv[y,s],
    C_liion = ω.values.C_liion[y,s],
    C_tes = ω.values.C_tes[y,s],
    C_tank = ω.values.C_tank[y,s],
    C_elyz = ω.values.C_elyz[y,s],
    C_fc = ω.values.C_fc[y,s],
    C_heater = ω.values.C_heater[y,s],
    # Electricity tariff
    C_grid_in = ω.values.C_grid_in[:,y,s],
    C_grid_out = ω.values.C_grid_out[:,y,s],
    )

    return Scenarios(ω.timestamp,values, [ω.proba[s]])
end
# Multistage designer reduction
function scenarios_reduction(designer::AbstractMultiStageDesigner, ω::Scenarios)

    if designer.options.scenario_reduction == "manual"
        # Scenario index is manually chosen
        s = designer.options.s
    else
        # Scenario index is randomly chosen to 1
        s = 1
    end

    # Reduced values
    values = (
    ld_E = ω.values.ld_E[:,:,s],
    ld_H = ω.values.ld_H[:,:,s],
    # Production
    pv_E = ω.values.pv_E[:,:,s],
    # Investment costs
    C_pv = ω.values.C_pv[:,s],
    C_liion = ω.values.C_liion[:,s],
    C_tes = ω.values.C_tes[:,s],
    C_tank = ω.values.C_tank[:,s],
    C_elyz = ω.values.C_elyz[:,s],
    C_fc = ω.values.C_fc[:,s],
    C_heater = ω.values.C_heater[:,s],
    # Electricity tariff
    C_grid_in = ω.values.C_grid_in[:,:,s],
    C_grid_out = ω.values.C_grid_out[:,:,s],
    )

    return Scenarios(ω.timestamp,values, [ω.proba[s]])
end
# Multistage stochastic designer reduction
function scenarios_reduction(designer::AbstractMultiStageStochasticDesigner, ω::Scenarios)

    if designer.options.scenario_reduction == "manual"
        # Scenario indexes are manually chosen
        s = designer.options.range_s
    else
        # TODO !! so far, scenario indexes are randomly chosen
        s = size(ω.values.ld_E,3)
    end

    # Reduced values
    values = (
    ld_E = ω.values.ld_E[:,:,s],
    ld_H = ω.values.ld_H[:,:,s],
    # Production
    pv_E = ω.values.pv_E[:,:,s],
    # Investment costs
    C_pv = ω.values.C_pv[:,s],
    C_liion = ω.values.C_liion[:,s],
    C_tes = ω.values.C_tes[:,s],
    C_tank = ω.values.C_tank[:,s],
    C_elyz = ω.values.C_elyz[:,s],
    C_fc = ω.values.C_fc[:,s],
    C_heater = ω.values.C_heater[:,s],
    # Electricity tariff
    C_grid_in = ω.values.C_grid_in[:,:,s],
    C_grid_out = ω.values.C_grid_out[:,:,s],
    )

    return Scenarios(ω.timestamp,values, [ω.proba[s]])
end
