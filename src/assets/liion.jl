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
     bounds::NamedTuple{(:lb, :ub), Tuple{Float64, Float64}}
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
        bounds = (lb = 0., ub = 1000.),
        Erated_ini = 1e-6,
        soc_ini = 0.5,
        soh_ini = 1.) =
        new(α_p_ch, α_p_dch, η_ch, η_dch, η_self, α_soc_min, α_soc_max, lifetime, nCycle, bounds, Erated_ini, soc_ini, soh_ini)
 end

### Preallocation
 function preallocate!(liion::Liion, nh::Int64, ny::Int64, ns::Int64)
     liion.Erated = convert(SharedArray,zeros(ny+1, ns)) ; liion.Erated[1,:] .= liion.Erated_ini
     liion.carrier = Electricity()
     liion.carrier.power = convert(SharedArray,zeros(nh, ny, ns))
     liion.soc = convert(SharedArray,zeros(nh+1, ny+1, ns)) ; liion.soc[1,1,:] .= liion.soc_ini
     liion.soh = convert(SharedArray,zeros(nh+1, ny+1, ns)) ; liion.soh[1,1,:] .= liion.soh_ini
     liion.cost = convert(SharedArray,zeros(ny, ns))
     return liion
 end

 ### Operation dynamic
function compute_operation_dynamics!(h::Int64, y::Int64, s::Int64, liion::Liion, decision::Float64, Δh::Int64)
     liion.soc[h+1,y,s], liion.soh[h+1,y,s], liion.carrier.power[h,y,s] = compute_operation_dynamics(liion, (Erated = liion.Erated[y,s], soc = liion.soc[h,y,s], soh = liion.soh[h,y,s]), decision, Δh)
end

function compute_operation_dynamics(liion::Liion, state::NamedTuple{(:Erated, :soc, :soh), Tuple{Float64, Float64, Float64}}, decision::Float64, Δh::Int64)
     # Control power constraint and correction
     power_dch = max(min(decision, liion.α_p_dch * state.Erated, state.soh * state.Erated / Δh, liion.η_dch * (state.soc * (1. - liion.η_self * Δh) - liion.α_soc_min) * state.Erated / Δh), 0.)
     power_ch = min(max(decision, -liion.α_p_ch * state.Erated, -state.soh * state.Erated / Δh, (state.soc * (1. - liion.η_self * Δh) - liion.α_soc_max) * state.Erated / Δh / liion.η_ch), 0.)
     # SoC dynamic
     soc_next = state.soc * (1. - liion.η_self * Δh) - (power_ch * liion.η_ch + power_dch / liion.η_dch) * Δh / state.Erated
     # SoH dynamic
     soh_next = state.soh - (power_dch - power_ch) * Δh / (2. * liion.nCycle * (liion.α_soc_max - liion.α_soc_min) * state.Erated)
     return soc_next, soh_next, power_dch + power_ch
end

 ### Investment dynamic
 function compute_investment_dynamics!(y::Int64, s::Int64, liion::Liion, decision::Union{Float64, Int64})
     liion.Erated[y+1,s], liion.soc[1,y+1,s], liion.soh[1,y+1,s] = compute_investment_dynamics(liion, (Erated = liion.Erated[y,s], soc = liion.soc[end,y,s], soh = liion.soh[end,y,s]), decision)
 end

 function compute_investment_dynamics(liion::Liion, state::NamedTuple{(:Erated, :soc, :soh), Tuple{Float64, Float64, Float64}}, decision::Union{Float64, Int64})
     if decision > 1e-2
         Erated_next = decision
         soc_next = liion.soc_ini
         soh_next =  1.
     else
         Erated_next = state.Erated
         soc_next = state.soc
         soh_next = state.soh
     end
     return Erated_next, soc_next, soh_next
 end
