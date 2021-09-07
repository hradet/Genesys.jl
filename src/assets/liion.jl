#=
    Li-ion battery modelling
 =#

 mutable struct Liion <: AbstractStorage
     # Parameters
     α_p_ch::Float64
     α_p_dch::Float64
     η_ch::Float64
     η_dch::Float64
     η_self::Float64
     α_soc_min::Float64
     α_soc_max::Float64
     lifetime::Int64
     nCycle::Float64
     # Initial conditions
     Erated_ini::Float64
     soc_ini::Float64
     soh_ini::Float64
     # Variables
     Erated::AbstractArray{Float64,2}
     carrier::Electricity
     soc::AbstractArray{Float64,3}
     soh::AbstractArray{Float64,3}
     # Eco
     cost::AbstractArray{Float64,2}
     # Inner constructor
     Liion(; α_p_ch = 1.5,
        α_p_dch = 1.5,
        η_ch = 0.9,
        η_dch = 0.9,
        η_self = 0.0005,
        α_soc_min = 0.2,
        α_soc_max = 0.8,
        lifetime = 12,
        nCycle = 2500.,
        Erated_ini = 1e-6,
        soc_ini = 0.5,
        soh_ini = 1.) =
        new(α_p_ch, α_p_dch, η_ch, η_dch, η_self, α_soc_min, α_soc_max, lifetime, nCycle, Erated_ini, soc_ini, soh_ini)
 end

### Preallocation
 function preallocate!(liion::Liion, nh::Int64, ny::Int64, ns::Int64)
     liion.Erated = convert(SharedArray,zeros(ny+1, ns)) ; liion.Erated[1,:] .= liion.Erated_ini
     liion.carrier = Electricity()
     liion.carrier.in = convert(SharedArray,zeros(nh, ny, ns))
     liion.carrier.out = convert(SharedArray,zeros(nh, ny, ns))
     liion.soc = convert(SharedArray,zeros(nh+1, ny+1, ns)) ; liion.soc[1,1,:] .= liion.soc_ini
     liion.soh = convert(SharedArray,zeros(nh+1, ny+1, ns)) ; liion.soh[1,1,:] .= liion.soh_ini
     liion.cost = convert(SharedArray,zeros(ny, ns))
     return liion
 end

 ### Operation dynamic
function compute_operation_dynamics!(h::Int64, y::Int64, s::Int64, liion::Liion, decision::Float64, Δh::Int64)
     # Control power constraint and correction
     liion.carrier.in[h,y,s] = max(min(decision, liion.α_p_dch * liion.Erated[y,s], liion.soh[h,y,s] * liion.Erated[y,s] / Δh, liion.η_dch * (liion.soc[h,y,s] * (1. - liion.η_self * Δh) - liion.α_soc_min) * liion.Erated[y,s] / Δh), 0.)
     liion.carrier.out[h,y,s] = min(max(decision, -liion.α_p_ch * liion.Erated[y,s], -liion.soh[h,y,s] * liion.Erated[y,s] / Δh, (liion.soc[h,y,s] * (1. - liion.η_self * Δh) - liion.α_soc_max) * liion.Erated[y,s] / Δh / liion.η_ch), 0.)
     # SoC dynamic
     liion.soc[h+1,y,s] = liion.soc[h,y,s] * (1. - liion.η_self * Δh) - (liion.carrier.out[h,y,s] * liion.η_ch + liion.carrier.in[h,y,s] / liion.η_dch) * Δh / liion.Erated[y,s]
     # SoH dynamic
     liion.soh[h+1,y,s] = liion.soh[h,y,s] - (liion.carrier.in[h,y,s] - liion.carrier.out[h,y,s]) * Δh / (2. * liion.nCycle * (liion.α_soc_max - liion.α_soc_min) * liion.Erated[y,s])
end

 ### Investment dynamic
 function compute_investment_dynamics!(y::Int64, s::Int64, liion::Liion, decision::Union{Float64, Int64})
     if decision > 1e-2
         liion.Erated[y+1,s] = decision
         liion.soc[1,y+1,s] = liion.soc[1,1,1]
         liion.soh[1,y+1,s] =  1.
     else
         liion.Erated[y+1,s] = liion.Erated[y,s]
         liion.soc[1,y+1,s] = liion.soc[end,y,s]
         liion.soh[1,y+1,s] = liion.soh[end,y,s]
     end
 end
