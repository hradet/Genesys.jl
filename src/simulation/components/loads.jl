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
function Load(outputGUI, nh, ny, ns)
     power_E = convert(SharedArray,zeros(nh, ny, ns))
     power_H = convert(SharedArray,zeros(nh, ny, ns))
     return Load(power_E,power_H)
end

 #                               Operation dynamic
 #______________________________________________________________________________


 #                               Investment dynamic
 #______________________________________________________________________________
