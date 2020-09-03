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
# Scenario Reduction function
function scenarios_reduction(ω::Scenarios, mode::String)

    if mode == "eac"
        ω_red = eac_reduction(ω)
    elseif mode == "eac_stoch"
        #TODO
    elseif mode == "metaheuristic"
        ω_red = metaheuristic_reduction(ω)
    elseif mode == "metaheuristic_stoch"
        #TODO
    end

    return ω_red
end
# Scenario reduction for EAC method
function eac_reduction(ω::Scenarios)
    # Only choose one year of data among the input scenarios
    y, s = 1, 1

    # Reduced values
    values = (
    ld_E = ω.values.ld_E[:,y,s],
    ld_H = ω.ld_H[:,y,s],
    # Production
    pv_E = ω.pv_E[:,y,s],
    # Investment costs
    C_pv = ω.C_pv[y,s],
    C_liion = ω.C_liion[y,s],
    C_tes = ω.C_tes[y,s],
    C_tank = ω.C_tank[y,s],
    C_elyz = ω.C_elyz[y,s],
    C_fc = ω.C_fc[y,s],
    C_heater = ω.C_heater[y,s],
    # Electricity tariff
    C_grid_in = ω.C_grid_in[:,y,s],
    C_grid_out = ω.C_grid_out[:,y,s],
    )

    return Scenarios(ω.timestamp,values, [ω.proba[s]])
end
# Scenario reduction for metaheuristic method
function metaheuristic_reduction(ω::Scenarios)
    # Only choose one scenario among the input scenarios
    s = 1

    # Reduced values
    values = (
    ld_E = ω.values.ld_E[:,:,s],
    ld_H = ω.ld_H[:,:,s],
    # Production
    pv_E = ω.pv_E[:,:,s],
    # Investment costs
    C_pv = ω.C_pv[:,s],
    C_liion = ω.C_liion[:,s],
    C_tes = ω.C_tes[:,s],
    C_tank = ω.C_tank[:,s],
    C_elyz = ω.C_elyz[:,s],
    C_fc = ω.C_fc[:,s],
    C_heater = ω.C_heater[:,s],
    # Electricity tariff
    C_grid_in = ω.C_grid_in[:,:,s],
    C_grid_out = ω.C_grid_out[:,:,s],
    )

    return Scenarios(ω.timestamp,values, [ω.proba[s]])
end
