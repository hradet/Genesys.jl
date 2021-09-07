#=
    Electrolyzer modelling
 =#

mutable struct Electrolyzer <: AbstractConverter
     # Paramètres
     α_p::Float64
     η_E_H2::Float64
     η_E_H::Float64
     lifetime::Int64
     nHoursMax::Float64
     # Initial conditions
     powerMax_ini::Float64
     soh_ini::Float64
     # Variables
     powerMax::AbstractArray{Float64,2}
     carrier::Vector{EnergyCarrier}
     soh::AbstractArray{Float64,3}
     # Eco
     cost::AbstractArray{Float64,2}
     # Inner constructor
     Electrolyzer(; α_p = 5/100,
                 η_E_H2 = 0.5,
                 η_E_H = 0.3,
                 lifetime = 15,
                 nHoursMax = 26000.,
                 powerMax_ini = 0.,
                 soh_ini = 1.) =
                 new(α_p, η_E_H2, η_E_H, lifetime, nHoursMax, powerMax_ini, soh_ini)
end

### Preallocation
function preallocate!(elyz::Electrolyzer, nh::Int64, ny::Int64, ns::Int64)
     elyz.powerMax = convert(SharedArray,zeros(ny+1, ns)) ; elyz.powerMax[1,:] .= elyz.powerMax_ini
     elyz.carrier = [Electricity(), Heat(), Hydrogen()]
     elyz.carrier[1].in = convert(SharedArray,zeros(nh, ny, ns))
     elyz.carrier[1].out = convert(SharedArray,zeros(nh, ny, ns))
     elyz.carrier[2].in = convert(SharedArray,zeros(nh, ny, ns))
     elyz.carrier[2].out = convert(SharedArray,zeros(nh, ny, ns))
     elyz.carrier[3].in = convert(SharedArray,zeros(nh, ny, ns))
     elyz.carrier[3].out = convert(SharedArray,zeros(nh, ny, ns))
     elyz.soh = convert(SharedArray,zeros(nh+1, ny+1, ns)) ; elyz.soh[1,1,:] .= elyz.soh_ini
     elyz.cost = convert(SharedArray,zeros(ny, ns))
     return elyz
end

### Operation dynamic
function compute_operation_dynamics!(h::Int64, y::Int64, s::Int64, elyz::Electrolyzer, decision::Float64, Δh::Int64)
    # Power constraint and correction
    elyz.α_p * elyz.powerMax[y,s] >= decision && elyz.soh[h,y,s] * elyz.nHoursMax / Δh > 1. ? elyz.carrier[1].out[h,y,s] = max(decision, -elyz.powerMax[y,s]) : elyz.carrier[1].out[h,y,s] = 0.
    # Power conversion
    elyz.carrier[3].in[h,y,s] = - elyz.carrier[1].out[h,y,s] * elyz.η_E_H2
    elyz.carrier[2].in[h,y,s] = - elyz.carrier[1].out[h,y,s] * elyz.η_E_H
    # SoH computation
    elyz.soh[h+1,y,s] = elyz.soh[h,y,s] - (elyz.carrier[1].out[h,y,s] > 0.) * Δh / elyz.nHoursMax
end

### Investment dynamic
function compute_investment_dynamics!(y::Int64, s::Int64, elyz::Electrolyzer, decision::Union{Float64, Int64})
    if decision > 1e-2
        elyz.powerMax[y+1,s] = decision
        elyz.soh[1,y+1,s] = 1.
    else
        elyz.powerMax[y+1,s] = elyz.powerMax[y,s]
        elyz.soh[1,y+1,s] = elyz.soh[end,y,s]
    end
end
