#=
    Electrolyzer modelling
 =#

mutable struct Electrolyzer
     # Paramètres
     α_p::Float64
     η_E_H2::Float64
     η_E_H::Float64
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
     C_elyz::AbstractArray{Float64,2}
     # Inner constructor
     Electrolyzer(; α_p = 5/100,
                 η_E_H2 = 0.5,
                 η_E_H = 0.3,
                 lifetime = 10.,
                 nHoursMax = 26000.,
                 powerMax_ini = 0.,
                 soh_ini = 1.) =
                 new(α_p, η_E_H2, η_E_H, lifetime, nHoursMax, powerMax_ini, soh_ini)
end

### Preallocation
function preallocate!(elyz::Electrolyzer, nh::Int64, ny::Int64, ns::Int64)
     elyz.powerMax = convert(SharedArray,zeros(ny+1, ns)) ; elyz.powerMax[1,:] .= elyz.powerMax_ini
     elyz.power_E = convert(SharedArray,zeros(nh, ny, ns))
     elyz.power_H = convert(SharedArray,zeros(nh, ny, ns))
     elyz.power_H2 = convert(SharedArray,zeros(nh, ny, ns))
     elyz.soh = convert(SharedArray,zeros(nh+1, ny+1, ns)) ; elyz.soh[1,1,:] .= elyz.soh_ini
     elyz.C_elyz = convert(SharedArray,zeros(ny, ns))
end

### Operation dynamic
function compute_operation_dynamics(elyz::Electrolyzer, x_elyz::NamedTuple, u_elyz::Float64, Δh::Int64)
    #=
    INPUT :
            x_elyz = (powerMax[y], soh[h,y]) tuple
            u_elyz[h,y] = control electric power in kW
    OUTPUT :
            power_E
            power_H
            power_H2
    =#

    # Power constraint and correction
    elyz.α_p * x_elyz.powerMax >= u_elyz && x_elyz.soh * elyz.nHoursMax / Δh > 1. ? power_E = max(u_elyz, -x_elyz.powerMax) : power_E = 0

    # Power computations
    power_H2 = - power_E * elyz.η_E_H2
    power_H = - power_E * elyz.η_E_H

    # SoH computation
    soh_next = x_elyz.soh - (power_E < 0.) * Δh / elyz.nHoursMax

    return power_E, power_H, power_H2, soh_next
end

### Investment dynamic
function compute_investment_dynamics(elyz::Electrolyzer, x_elyz::NamedTuple, u_elyz::Union{Float64, Int64})
    #=
        INPUT :
                x_elyz = [powerMax[y], soh[end,y]]
                u_elyz[y] = elyz control inv in kW
        OUTPUT :
                pMax_next
                soh_next
    =#

    # Model
    if u_elyz > 1e-2
        powerMax_next = u_elyz
        soh_next = 1.
    else
        powerMax_next = x_elyz.powerMax
        soh_next = x_elyz.soh
    end

    return powerMax_next, soh_next
end
