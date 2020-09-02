#=
    Sources modelling
 =#

 #                                  Structure
 #______________________________________________________________________________
mutable struct Source
     lifetime::Float64
     # Variables
     power_E::AbstractArray{Float64,3}
     powerMax::AbstractArray{Float64,2}
     # Eco
     C_pv::AbstractArray{Float64,2}
end
# Constructor
function Source(outputGUI::NamedTuple, nh::Int64, ny::Int64, ns::Int64)
     lifetime = outputGUI.lifetime
     power_E = convert(SharedArray,zeros(nh, ny, ns))
     powerMax = convert(SharedArray,zeros(ny+1, ns))
     C_pv = convert(SharedArray,zeros(ny, ns))
     return Source(lifetime,power_E,powerMax,C_pv)
 end

 #                               Operation dynamic
 #______________________________________________________________________________


 #                               Investment dynamic
 #______________________________________________________________________________
 function compute_investment_dynamics(pv::Source, x_pv::NamedTuple, u_pv::Float64)
     #=
         INPUT :
                 x_pv = [powerMax[y]]
                 u_pv[y] = pv control inv in kW
         OUTPUT :
                 pMax_next
     =#

     # Model
     if round(u_pv) > 0.
         powerMax_next = u_pv
     else
         powerMax_next = x_pv.powerMax
     end

     return powerMax_next
 end
