#=
    Heater modelling
 =#

mutable struct Heater <: AbstractConverter
     # Parameters
     η_E_H::Float64
     lifetime::Int64
     # Initial conditions
     powerMax_ini::Float64
     soh_ini::Float64
     # Variables
     powerMax::AbstractArray{Float64,2}
     carrier::Vector{Any}
     # Eco
     cost::AbstractArray{Float64,2}
     # Inner constructor
     Heater(; η_E_H = 1.,
             lifetime = 25.,
             powerMax_ini = 30.,
             soh_ini = 1.) =
             new(η_E_H, lifetime, powerMax_ini)
end

### Preallocation
function preallocate!(heater::Heater, nh::Int64, ny::Int64, ns::Int64)
     heater.powerMax = convert(SharedArray,zeros(ny+1, ns)) ; heater.powerMax[1,:] .= heater.powerMax_ini
     heater.carrier = [Electricity(), Heat()]
     heater.carrier[1].in = convert(SharedArray,zeros(nh, ny, ns))
     heater.carrier[1].out = convert(SharedArray,zeros(nh, ny, ns))
     heater.carrier[2].in = convert(SharedArray,zeros(nh, ny, ns))
     heater.carrier[2].out = convert(SharedArray,zeros(nh, ny, ns))
     heater.cost = convert(SharedArray,zeros(ny, ns))
     return heater
end

 ### Operation dynamic
 function compute_operation_dynamics!(h::Int64, y::Int64, s::Int64, heater::Heater, decision::Float64, Δh::Int64)
     # Power constraint and correction
     heater.carrier[1].out[h,y,s] = min(max(decision, -heater.powerMax[y,s]), 0.)
     # Power computation
     heater.carrier[2].in[h,y,s] = - heater.η_E_H * heater.carrier[1].out[h,y,s]
 end

 ### Investment dynamic
 function compute_investment_dynamics!(y::Int64, s::Int64, heater::Heater, decision::Union{Float64, Int64})
     if decision > 1e-2
         heater.powerMax[y+1,s] = decision
     else
         heater.powerMax[y+1,s] = heater.powerMax[y,s]
     end
 end
