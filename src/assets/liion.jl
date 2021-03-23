#=
    Li-ion battery modelling
 =#

 mutable struct Liion
     # Parameters
     α_p_ch::Float64
     α_p_dch::Float64
     η_ch::Float64
     η_dch::Float64
     η_self::Float64
     α_soc_min::Float64
     α_soc_max::Float64
     lifetime::Float64
     nCycle::Float64
     # Initial conditions
     Erated_ini::Float64
     soc_ini::Float64
     soh_ini::Float64
     # Variables
     Erated::AbstractArray{Float64,2}
     power_E::AbstractArray{Float64,3}
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
        nCycle = 2500,
        Erated_ini = 1e-6,
        soc_ini = 0.5,
        soh_ini = 1.) =
        new(α_p_ch, α_p_dch, η_ch, η_dch, η_self, α_soc_min, α_soc_max, lifetime, nCycle, Erated_ini, soc_ini, soh_ini)
 end

### Preallocation
 function preallocate!(liion::Liion, nh::Int64, ny::Int64, ns::Int64)
     liion.Erated = convert(SharedArray,zeros(ny+1, ns)) ; liion.Erated[1,:] .= liion.Erated_ini
     liion.power_E = convert(SharedArray,zeros(nh, ny, ns))
     liion.soc = convert(SharedArray,zeros(nh+1, ny+1, ns)) ; liion.soc[1,1,:] .= liion.soc_ini
     liion.soh = convert(SharedArray,zeros(nh+1, ny+1, ns)) ; liion.soh[1,1,:] .= liion.soh_ini
     liion.cost = convert(SharedArray,zeros(ny, ns))
 end

 ### Operation dynamic
function compute_operation_dynamics(liion::Liion, x_liion::NamedTuple{(:Erated, :soc, :soh), Tuple{Float64, Float64, Float64}}, u_liion::Float64, Δh::Int64)
     #=
     INPUT :
             x_liion = (Erated[y], soc[h,y], soh[h,y]) tuple
             u_liion[h,y] = control power in kW
     OUTPUT :
             soc_next
             soh_next
             power = the real battery power in kW
     =#

     # Control power constraint and correction
     power_dch = max(min(u_liion, liion.α_p_dch * x_liion.Erated, x_liion.soh * x_liion.Erated / Δh, liion.η_dch * (x_liion.soc * (1. - liion.η_self * Δh) - liion.α_soc_min) * x_liion.Erated / Δh), 0.)
     power_ch = min(max(u_liion, -liion.α_p_ch * x_liion.Erated, -x_liion.soh * x_liion.Erated / Δh, (x_liion.soc * (1. - liion.η_self * Δh) - liion.α_soc_max) * x_liion.Erated / Δh / liion.η_ch), 0.)

     # SoC dynamic
     soc_next = x_liion.soc * (1. - liion.η_self * Δh) - (power_ch * liion.η_ch + power_dch / liion.η_dch) * Δh / x_liion.Erated

     # SoH dynamic
     soh_next = x_liion.soh - (power_dch - power_ch) * Δh / (2. * liion.nCycle * (liion.α_soc_max - liion.α_soc_min) * x_liion.Erated)

     return soc_next, soh_next, power_ch + power_dch
end

 ### Investment dynamic
 function compute_investment_dynamics(liion::Liion, x_liion::NamedTuple{(:Erated, :soc, :soh), Tuple{Float64, Float64, Float64}}, u_liion::Union{Float64, Int64})
     #=
         INPUT :
                 x_liion = [Erated[y], soc[end,y], soh[end,y]]
                 u_liion[y] = liion control inv in kWh
         OUTPUT :
                 E_next
                 soc_next
                 soh_next
     =#

     # Model
     if u_liion > 1e-2
         E_next = u_liion
         soc_next = liion.soc[1,1,1]
         soh_next =  1.
     else
         E_next = x_liion.Erated
         soc_next = x_liion.soc
         soh_next =  x_liion.soh
     end

     return E_next, soc_next, soh_next
 end
