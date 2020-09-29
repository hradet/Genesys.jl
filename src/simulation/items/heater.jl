#=
    Heater modelling
 =#

mutable struct Heater
     # Parameters
     η_E_H::Float64
     lifetime::Float64
     # Initial conditions
     powerMax_ini
     # Variables
     powerMax::AbstractArray{Float64,2}
     power_E::AbstractArray{Float64,3}
     power_H::AbstractArray{Float64,3}
     # Eco
     C_heater::AbstractArray{Float64,2}
     # Inner constructor
     Heater(; η_E_H = 1.,
             lifetime = 25.,
             powerMax_ini = 30.) =
             new(η_E_H, lifetime, powerMax_ini)
end

### Preallocation
function preallocate!(heater::Heater, nh::Int64, ny::Int64, ns::Int64)
     heater.powerMax = convert(SharedArray,zeros(ny+1, ns)) ; heater.powerMax .= heater.powerMax_ini
     heater.power_E = convert(SharedArray,zeros(nh, ny, ns))
     heater.power_H = convert(SharedArray,zeros(nh, ny, ns))
     heater.C_heater = convert(SharedArray,zeros(ny, ns))
end

 ### Operation dynamic
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
     0. <= -u_heater <= x_heater.powerMax ? power_E = u_heater : power_E = 0.

     # Power computation
     power_H = - heater.η_E_H * power_E

     return power_E, power_H
 end

 ### Investment dynamic
