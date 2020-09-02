#=
    Heater modelling
 =#

#                                 Structure
#______________________________________________________________________________
struct Heater
     # Parameters
     η_E_H::Float64
     lifetime::Float64
     # Variables
     powerMax::AbstractArray{Float64,2}
     power_E::AbstractArray{Float64,3}
     power_H::AbstractArray{Float64,3}
     # Eco
     C_heater::AbstractArray{Float64,2}
end
# Constructor
function Heater(outputGUI::NamedTuple, nh::Int64, ny::Int64, ns::Int64)
     # Parameters
     η_E_H = outputGUI.η_E_H
     lifetime = outputGUI.lifetime
     # Variables
     powerMax = convert(SharedArray,zeros(ny+1, ns))
     power_E = convert(SharedArray,zeros(nh, ny, ns))
     power_H = convert(SharedArray,zeros(nh, ny, ns))
     # Eco
     C_heater = convert(SharedArray,zeros(ny, ns))
     return Heater(η_E_H,lifetime,powerMax,power_E,power_H,C_heater)
end

 #                               Operation dynamic
 #______________________________________________________________________________
 function compute_operation_dynamics(heater::Heater, x_heater::NamedTuple, u_heater::Float64, Δh::Int64)
     #=
     INPUT :
             x_heater = (powerMax[y]) tuple
             u_heater[h,y] = control electric power in kW
     OUTPUT :
             power_E
             power_H
     =#

     # Power constraint and correction
     0. <= -u_heater <= x_heater.powerMax ? power_E = u_heater : power_E = 0

     # Power computation
     power_H = - heater.η_E_H * power_E

     return power_E, power_H
 end

 #                               Investment dynamic
 #______________________________________________________________________________
