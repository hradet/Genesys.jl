#=
    Fuel cell modelling
 =#

mutable struct FuelCell <: AbstractConverter
     # Paramètres
     α_p::Float64
     η_H2_E::Float64
     η_H2_H::Float64
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
     FuelCell(; α_p = 8/100,
             η_H2_E = 0.4,
             η_H2_H = 0.4,
             lifetime = 14,
             nHoursMax = 10000.,
             powerMax_ini = 0.,
             soh_ini = 1.) =
             new(α_p, η_H2_E, η_H2_H, lifetime, nHoursMax, powerMax_ini, soh_ini)
end

### Preallocation
function preallocate!(fc::FuelCell, nh::Int64, ny::Int64, ns::Int64)
     fc.powerMax = convert(SharedArray,zeros(ny+1, ns)) ; fc.powerMax[1,:] .= fc.powerMax_ini
     fc.carrier = [Electricity(), Heat(), Hydrogen()]
     fc.carrier[1].in = convert(SharedArray,zeros(nh, ny, ns))
     fc.carrier[1].out = convert(SharedArray,zeros(nh, ny, ns))
     fc.carrier[2].in = convert(SharedArray,zeros(nh, ny, ns))
     fc.carrier[2].out = convert(SharedArray,zeros(nh, ny, ns))
     fc.carrier[3].in = convert(SharedArray,zeros(nh, ny, ns))
     fc.carrier[3].out = convert(SharedArray,zeros(nh, ny, ns))
     fc.soh = convert(SharedArray,zeros(nh+1, ny+1, ns)) ; fc.soh[1,1,:] .= fc.soh_ini
     fc.cost = convert(SharedArray,zeros(ny, ns))
     return fc
end

### Operation dynamic
function compute_operation_dynamics!(h::Int64, y::Int64, s::Int64, fc::FuelCell, decision::Float64, Δh::Int64)
    # Power constraint and correction
    fc.α_p * fc.powerMax[y,s] <= decision && fc.soh[h,y,s] * fc.nHoursMax / Δh > 1. ? fc.carrier[1].in[h,y,s] = min(decision, fc.powerMax[y,s]) : fc.carrier[1].in[h,y,s] = 0.
    # Power conversion
    fc.carrier[3].out[h,y,s] = - fc.carrier[1].in[h,y,s] / fc.η_H2_E
    fc.carrier[2].in[h,y,s] = - fc.carrier[3].out[h,y,s] * fc.η_H2_H
    # SoH computation
    fc.soh[h+1,y,s] = fc.soh[h,y,s] - (fc.carrier[1].in[h,y,s] > 0.) * Δh / fc.nHoursMax
end

### Investment dynamic
function compute_investment_dynamics!(y::Int64, s::Int64, fc::FuelCell, decision::Union{Float64, Int64})
    if decision > 1e-2
        fc.powerMax[y+1,s] = decision
        fc.soh[1,y+1,s] = 1.
    else
        fc.powerMax[y+1,s] = fc.powerMax[y,s]
        fc.soh[1,y+1,s] = fc.soh[end,y,s]
    end
end
