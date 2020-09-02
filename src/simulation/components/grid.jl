#=
    Elec grid modelling
 =#

#                                  Structure
#______________________________________________________________________________
# TODO : Mettre les paramètres tau dans parameter + modifier optimization function
mutable struct Grid
     # Parameters
     τ_power
     τ_energy
     # Variables
     power_E::AbstractArray{Float64,3}
     C_grid_in::AbstractArray{Float64,3}
     C_grid_out::AbstractArray{Float64,3}
end
# Constructor
function Grid(outputGUI::NamedTuple, nh::Int64, ny::Int64, ns::Int64)
     # Parameters
     τ_power = outputGUI.τ_power
     τ_energy = outputGUI.τ_energy
     # Variables
     power_E = convert(SharedArray,zeros(nh, ny, ns))
     # Eco
     C_grid_in = convert(SharedArray,zeros(nh, ny, ns))
     C_grid_out = convert(SharedArray,zeros(nh, ny, ns))
     return Grid(τ_power,τ_energy,power_E,C_grid_in,C_grid_out)
end

#                               Operation dynamic
#______________________________________________________________________________


#                               Investment dynamic
#______________________________________________________________________________
