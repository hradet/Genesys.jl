#=
    Thermal energy storage modelling
 =#

mutable struct ThermalSto
     # Paramètres
     α_p_ch::Float64
     α_p_dch::Float64
     η_ch::Float64
     η_dch::Float64
     η_self::Float64
     α_soc_min::Float64
     α_soc_max::Float64
     lifetime::Float64
     # Initial conditions
     Erated_ini::Float64
     soc_ini::Float64
     # Variables
     Erated::AbstractArray{Float64,2}
     power_H::AbstractArray{Float64,3}
     soc::AbstractArray{Float64,3}
     # Eco
     C_tes::AbstractArray{Float64,2}
     # Inner constructor
     ThermalSto(; α_p_ch = 1.5,
          α_p_dch = 1.5,
          η_ch = 0.9,
          η_dch = 0.9,
          η_self = 0.01,
          α_soc_min = 0.,
          α_soc_max = 1.,
          lifetime = 25,
          Erated_ini = 0.,
          soc_ini = 0.5) =
          new(α_p_ch, α_p_dch, η_ch, η_dch, η_self, α_soc_min, α_soc_max, lifetime, Erated_ini, soc_ini)
end

### Preallocation
function preallocate!(tes::ThermalSto, nh::Int64, ny::Int64, ns::Int64)
   tes.Erated = convert(SharedArray,zeros(ny+1, ns)) ; tes.Erated[1,:] .= tes.Erated_ini
   tes.power_H = convert(SharedArray,zeros(nh, ny, ns))
   tes.soc = convert(SharedArray,zeros(nh+1, ny+1, ns)) ; tes.soc[1,1,:] .= tes.soc_ini
   tes.C_tes = convert(SharedArray,zeros(ny, ns))
end

### Operation dynamic
function compute_operation_dynamics(tes::ThermalSto, x_tes::NamedTuple, u_tes::Float64, Δh::Int64)
    #=
    INPUT :
            x_tes = (Erated[y], soc[h,y]) tuple
            u_tes[h,y] = control power in kW
    OUTPUT :
            soc_next
            power = the real battery power in kW
    =#

    # Power constraint and correction
    0. <= u_tes <= tes.α_p_dch * x_tes.Erated ? power_dch = u_tes : power_dch = 0.
    0. <= -u_tes <= tes.α_p_ch * x_tes.Erated ? power_ch = u_tes : power_ch = 0.

    # SoC dynamic
    soc_next = x_tes.soc * (1. - tes.η_self * Δh) - (power_ch * tes.η_ch + power_dch / tes.η_dch) * Δh / x_tes.Erated

    # State variable bounds
    overshoot = (round(soc_next;digits=3) < tes.α_soc_min) || (round(soc_next;digits=3) > tes.α_soc_max)

    overshoot ? soc_next = max(x_tes.soc * (1. - tes.η_self * Δh), tes.α_soc_min) : nothing
    overshoot ? power_ch = power_dch = 0. : nothing

    return soc_next, power_ch + power_dch
end

### Investment dynamic
function compute_investment_dynamics(tes::ThermalSto, x_tes::NamedTuple, u_tes::Union{Float64, Int64})
    #=
        INPUT :
                x_tes = [Erated[y], soc[end,y]]
                u_tes[y] = tes control inv in kWh
        OUTPUT :
                E_next
                soc_next
    =#

    # Model
    if round(u_tes) > 0.
        E_next = u_tes
        soc_next = tes.soc[1,1,1]
    else
        E_next = x_tes.Erated
        soc_next = x_tes.soc
    end

    return E_next, soc_next
end
