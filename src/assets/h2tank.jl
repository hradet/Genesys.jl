#=
    H2 tank storage modelling
 =#

mutable struct H2Tank
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
     # Variable
     Erated::AbstractArray{Float64,2}
     power_H2::AbstractArray{Float64,3}
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
        soc_ini = 0.5) =
        new(α_p_ch, α_p_dch, η_ch, η_dch, η_self, α_soc_min, α_soc_max, lifetime, Erated_ini, soc_ini)
end

### Preallocation
function preallocate!(h2tank::H2Tank, nh::Int64, ny::Int64, ns::Int64)
   h2tank.Erated = convert(SharedArray,zeros(ny+1, ns)) ; h2tank.Erated[1,:] .= h2tank.Erated_ini
   h2tank.power_H2 = convert(SharedArray,zeros(nh, ny, ns))
   h2tank.soc = convert(SharedArray,zeros(nh+1, ny+1, ns)) ; h2tank.soc[1,1,:] .= h2tank.soc_ini
   h2tank.cost = convert(SharedArray,zeros(ny, ns))
end

### Operation dynamic
function compute_operation_dynamics(h2tank::H2Tank, x_h2tank::NamedTuple{(:Erated, :soc), Tuple{Float64, Float64}}, u_h2tank::Float64, Δh::Int64)
     #=
     INPUT :
             x_h2 = (Erated[y], soc[h,y]) tuple
             u_h2[h,y] = control power in kW
     OUTPUT :
             soc_next
             power = the real battery power in kW
     =#

     # Power constraint and correction
     # Control power constraint and correction
      power_dch = max(min(u_h2tank, h2tank.α_p_dch * x_h2tank.Erated, h2tank.η_dch * (x_h2tank.soc * (1. - h2tank.η_self * Δh) - h2tank.α_soc_min) * x_h2tank.Erated / Δh), 0.)
      power_ch = min(max(u_h2tank, -h2tank.α_p_ch * x_h2tank.Erated, (x_h2tank.soc * (1. - h2tank.η_self * Δh) - h2tank.α_soc_max) * x_h2tank.Erated / Δh / h2tank.η_ch), 0.)

      # SoC dynamic
      soc_next = x_h2tank.soc * (1. - h2tank.η_self * Δh) - (power_ch * h2tank.η_ch + power_dch / h2tank.η_dch) * Δh / x_h2tank.Erated

     return soc_next, power_ch + power_dch
end

### Investment dynamic
function compute_investment_dynamics(h2tank::H2Tank, x_tank::NamedTuple{(:Erated, :soc), Tuple{Float64, Float64}}, u_tank::Union{Float64, Int64})
     #=
         INPUT :
                 x_tank = [Erated[y], soc[end,y]]
                 u_tank[y] = h2tank control inv in kWh
         OUTPUT :
                 E_next
                 soc_next
     =#

     # Model
     if u_tank > 1e-2
         E_next = u_tank
         soc_next = h2tank.soc[1,1,1]
     else
         E_next = x_tank.Erated
         soc_next = x_tank.soc
     end

     return E_next, soc_next
end
