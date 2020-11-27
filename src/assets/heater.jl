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
     cost::AbstractArray{Float64,2}
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
     heater.cost = convert(SharedArray,zeros(ny, ns))
end

 ### Operation dynamic
 function compute_operation_dynamics(heater::Heater, x_heater::NamedTuple{(:powerMax,), Tuple{Float64}}, u_heater::Float64, Δh::Int64)
     #=
     INPUT :
             x_heater = (powerMax[y]) tuple
             u_heater[h,y] = control electric power in kW
     OUTPUT :
             power_E
             power_H
     =#

     # Power constraint and correction
     power_E = min(max(u_heater, -x_heater.powerMax), 0.)

     # Power computation
     power_H = - heater.η_E_H * power_E

     return power_E, power_H
 end

 ### Investment dynamic
