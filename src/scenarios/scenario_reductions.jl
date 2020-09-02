#=
    Scenario reduction functions
=#


mutable struct Scenarios <: AbstractScenarios
     # Timestamp
     timestamp::Union{Array{Float64,1}, Array{DateTime,1}}
     # Demand
     ld_E::Array{Float64,3}
     ld_H::Array{Float64,3}
     # Production
     pv_E::Array{Float64,3}
     # Investment costs
     C_pv::Array{Float64,2}
     C_liion::Array{Float64,2}
     C_tes::Array{Float64,2}
     C_tank::Array{Float64,2}
     C_elyz::Array{Float64,2}
     C_fc::Array{Float64,2}
     C_heater::Array{Float64,2}
     # Electricity tariff
     C_grid_in::Array{Float64,3}
     C_grid_out::Array{Float64,3}
end

# Constructor
function Scenarios(outputGUI::Dict{}, nh::Int64, ny::Int64, ns::Int64)
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

function scenarios_reduction(ω::Scenarios; mode = "reference")

    if mode == "oneyear"

    elseif mode == "typicaldays"

    elseif mode == "multiscenario"

    elseif mode == "multiscenario_typicaldays"

    end


    return
end
