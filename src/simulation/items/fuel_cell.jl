#=
    Fuel cell modelling
 =#

mutable struct FuelCell
     # Paramètres
     α_p::Float64
     η_H2_E::Float64
     η_H2_H::Float64
     lifetime::Float64
     nHoursMax::Float64
     # Initial conditions
     powerMax_ini::Float64
     soh_ini::Float64
     # Variables
     powerMax::AbstractArray{Float64,2}
     power_E::AbstractArray{Float64,3}
     power_H::AbstractArray{Float64,3}
     power_H2::AbstractArray{Float64,3}
     soh::AbstractArray{Float64,3}
     # Eco
     C_fc::AbstractArray{Float64,2}
     # Inner constructor
     FuelCell(; α_p = 8/100,
             η_H2_E = 0.4,
             η_H2_H = 0.4,
             lifetime = 5.,
             nHoursMax = 10000.,
             powerMax_ini = 0.,
             soh_ini = 1.) =
             new(α_p, η_H2_E, η_H2_H, lifetime, nHoursMax, powerMax_ini, soh_ini)
end

### Preallocation
function preallocate!(fc::FuelCell, nh::Int64, ny::Int64, ns::Int64)
     fc.powerMax = convert(SharedArray,zeros(ny+1, ns)) ; fc.powerMax[1,:] .= fc.powerMax_ini
     fc.power_E = convert(SharedArray,zeros(nh, ny, ns))
     fc.power_H = convert(SharedArray,zeros(nh, ny, ns))
     fc.power_H2 = convert(SharedArray,zeros(nh, ny, ns))
     fc.soh = convert(SharedArray,zeros(nh+1, ny+1, ns))
     fc.C_fc = convert(SharedArray,zeros(ny, ns))
end

### Operation dynamic
function compute_operation_dynamics(fc::FuelCell, x_fc::NamedTuple, u_fc::Float64, Δh::Int64)
     #=
     INPUT :
             x_fc = (powerMax[y], soh[h,y]) tuple
             u_fc[h,y] = control electric power in kW
     OUTPUT :
             power_E
             power_H
             power_H2
     =#

     # Power constraint and correction
     fc.α_p * x_fc.powerMax <= u_fc <= x_fc.powerMax ? power_E = u_fc : power_E = 0.

     # Power computations
     power_H2 = - power_E / fc.η_H2_E
     power_H = - power_H2 * fc.η_H2_H

     # SoH computation
     soh_next = x_fc.soh - (power_E > 0.) * Δh / fc.nHoursMax

     # SoH constraint and correction
     soh_next < 0. ? power_E = power_H = power_H2 = 0. : nothing
     soh_next < 0. ? soh_next = x_fc.soh : nothing

     return power_E, power_H, power_H2, soh_next
end

### Investment dynamic
function compute_investment_dynamics(fc::FuelCell, x_fc::NamedTuple, u_fc::Union{Float64, Int64})
     #=
         INPUT :
                 x_fc = [powerMax[y], soh[end,y]]
                 u_fc[y] = fc control inv in kW
         OUTPUT :
                 pMax_next
                 soh_next
     =#

     # Model
     if round(u_fc) > 0.
         powerMax_next = u_fc
         soh_next = 1.
     else
         powerMax_next = x_fc.powerMax
         soh_next = x_fc.soh
     end

     return powerMax_next, soh_next
end
