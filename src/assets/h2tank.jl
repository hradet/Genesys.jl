#=
    H2 tank storage modelling
 =#

mutable struct H2Tank  <: AbstractStorage
     # Paramètres
     α_p_ch::Float64
     α_p_dch::Float64
     η_ch::Float64
     η_dch::Float64
     η_self::Float64
     α_soc_min::Float64
     α_soc_max::Float64
     lifetime::Int64
     # Initial conditions
     Erated_ini::Float64
     soc_ini::Float64
     soh_ini::Float64
     # Variable
     Erated::AbstractArray{Float64,2}
     carrier::Hydrogen
     soc::AbstractArray{Float64,3}
     # Eco
     cost::AbstractArray{Float64,2}
     # Inner constructor
     H2Tank(; α_p_ch = 1.5,
        α_p_dch = 1.5,
        η_ch = 1.,
        η_dch = 1.,
        η_self = 0.,
        α_soc_min = 0.,
        α_soc_max = 1.,
        lifetime = 25,
        Erated_ini = 1e-6,
        soc_ini = 0.5,
        soh_ini = 1.) =
        new(α_p_ch, α_p_dch, η_ch, η_dch, η_self, α_soc_min, α_soc_max, lifetime, Erated_ini, soc_ini, soh_ini)
end

### Preallocation
function preallocate!(h2tank::H2Tank, nh::Int64, ny::Int64, ns::Int64)
    h2tank.Erated = convert(SharedArray,zeros(ny+1, ns)) ; h2tank.Erated[1,:] .= h2tank.Erated_ini
    h2tank.carrier = Hydrogen()
    h2tank.carrier.in = convert(SharedArray,zeros(nh, ny, ns))
    h2tank.carrier.out = convert(SharedArray,zeros(nh, ny, ns))
    h2tank.soc = convert(SharedArray,zeros(nh+1, ny+1, ns)) ; h2tank.soc[1,1,:] .= h2tank.soc_ini
    h2tank.cost = convert(SharedArray,zeros(ny, ns))
    return h2tank
end

### Operation dynamic
function compute_operation_dynamics!(h::Int64, y::Int64, s::Int64, h2tank::H2Tank, decision::Float64, Δh::Int64)
     # Control power constraint and correction
     power_in = max(min(decision, h2tank.α_p_dch * h2tank.Erated[y,s], h2tank.η_dch * (h2tank.soc[h,y,s] * (1. - h2tank.η_self * Δh) - h2tank.α_soc_min) * h2tank.Erated[y,s] / Δh), 0.)
     power_out = min(max(decision, -h2tank.α_p_ch * h2tank.Erated[y,s], (h2tank.soc[h,y,s] * (1. - h2tank.η_self * Δh) - h2tank.α_soc_max) * h2tank.Erated[y,s] / Δh / h2tank.η_ch), 0.)
     # SoC dynamic
     h2tank.soc[h+1,y,s] = h2tank.soc[h,y,s] * (1. - h2tank.η_self * Δh) - (power_out * h2tank.η_ch + power_in / h2tank.η_dch) * Δh / h2tank.Erated[y,s]
end

### Investment dynamic
function compute_investment_dynamics!(y::Int64, s::Int64, h2tank::H2Tank, decision::Union{Float64, Int64})
    if decision > 1e-2
        h2tank.Erated[y+1,s] = decision
        h2tank.soc[1,y+1,s] = h2tank.soc[1,1,1]
    else
        h2tank.Erated[y+1,s] = h2tank.Erated[y,s]
        h2tank.soc[1,y+1,s] = h2tank.soc[end,y,s]
    end
end
