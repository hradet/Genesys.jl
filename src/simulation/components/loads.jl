#=
    Loads modelling
 =#

 #                                  Structure
 #______________________________________________________________________________
mutable struct Load
     # Variables
     power_E::AbstractArray{Float64,3}
     power_H::AbstractArray{Float64,3}
end
# Constructor
function Load(outputGUI::Dict{String,Any}, nh::Int64, ny::Int64, ns::Int64)
     power_E = convert(SharedArray,zeros(nh, ny, ns))
     power_H = convert(SharedArray,zeros(nh, ny, ns))
     return Load(power_E,power_H)
end

 #                               Operation dynamic
 #______________________________________________________________________________


 #                               Investment dynamic
 #______________________________________________________________________________
