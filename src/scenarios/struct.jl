mutable struct Scenario <: AbstractScenarios
    data
    stats
end

mutable struct Scenarios <: AbstractScenarios
     # Timestamp
     timestamp
     # Demand
     ld_E
     ld_H
     # Production
     pv_E
     # Investment costs
     C_pv
     C_liion
     C_tes
     C_tank
     C_elyz
     C_fc
     C_heater
     # Electricity tariff
     C_grid_in
     C_grid_out
end
# Constructor
function Scenarios(outputGUI, nh, ny, ns)
     # Demand
     ld_E = outputGUI["ld_E"][1:nh, 1:ny, 1:ns]
     ld_H =outputGUI["ld_H"][1:nh, 1:ny, 1:ns]
     # Production
     pv_E = outputGUI["pv_E"][1:nh, 1:ny, 1:ns]
     # Investment costs
     C_pv = outputGUI["C_pv"][1:ny, 1:ns]
     C_liion = outputGUI["C_liion"][1:ny, 1:ns]
     C_tes = outputGUI["C_tes"][1:ny, 1:ns]
     C_tank = outputGUI["C_tank"][1:ny, 1:ns]
     C_elyz = outputGUI["C_elyz"][1:ny, 1:ns]
     C_fc = outputGUI["C_fc"][1:ny, 1:ns]
     C_heater = outputGUI["C_heater"][1:ny, 1:ns]
     # Electricity tariff
     C_grid_in = outputGUI["C_grid_in"][1:nh, 1:ny, 1:ns]
     C_grid_out = outputGUI["C_grid_out"][1:nh, 1:ny, 1:ns]
     return Scenarios(outputGUI["timestamp"],ld_E,ld_H,pv_E,C_pv,C_liion,C_tes,C_tank,C_elyz,C_fc,C_heater,C_grid_in,C_grid_out)
end

mutable struct ClusteredScenarios
    clusters::Scenarios # clusters
    σ # sequence of clusters = assignments
    nby # number of entities by cluster = counts
end

mutable struct MarkovChain
    states
    transition_matrices
end
