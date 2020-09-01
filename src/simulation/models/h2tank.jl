#=
    H2 tank storage modelling
 =#

 #                                  Structure
 #______________________________________________________________________________
struct H2Tank
     # Paramètres
     α_p_ch::Float64
     α_p_dch::Float64
     η_ch::Float64
     η_dch::Float64
     η_self::Float64
     α_soc_min::Float64
     α_soc_max::Float64
     lifetime::Float64
     # Variable
     Erated::AbstractArray{Float64,2}
     power_H2::AbstractArray{Float64,3}
     soc::AbstractArray{Float64,3}
     # Eco
     C_tank::AbstractArray{Float64,2}
end
# Constructor
function H2Tank(outputGUI, nh, ny, ns)
     # Paramètres
     α_p_ch = outputGUI.α_p_ch
     α_p_dch = outputGUI.α_p_dch
     η_ch = outputGUI.η_ch
     η_dch = outputGUI.η_dch
     η_self = outputGUI.η_self
     α_soc_min = outputGUI.α_soc_min
     α_soc_max = outputGUI.α_soc_max
     lifetime = outputGUI.lifetime
     # Variables
     Erated = convert(SharedArray,zeros(ny+1, ns))
     power_H2 = convert(SharedArray,zeros(nh, ny, ns))
     soc = convert(SharedArray,zeros(nh+1, ny+1, ns))
     # Eco
     C_tank = convert(SharedArray,zeros(ny, ns))
     return H2Tank(α_p_ch,α_p_dch,η_ch,η_dch,η_self,α_soc_min,α_soc_max,lifetime,Erated,power_H2,soc,C_tank)
end

#                               Operation dynamic
#______________________________________________________________________________
function compute_operation_dynamics(h2tank::H2Tank, x_h2::NamedTuple, u_h2::Float64, Δh::Int64)
     #=
     INPUT :
             x_h2 = (Erated[y], soc[h,y]) tuple
             u_h2[h,y] = control power in kW
     OUTPUT :
             soc_next
             power = the real battery power in kW
     =#

     # Power constraint and correction
     0 <= u_h2 <= h2tank.α_p_dch * x_h2.Erated ? power_dch = u_h2 : power_dch = 0
     0 <= -u_h2 <= h2tank.α_p_ch * x_h2.Erated ? power_ch = u_h2 : power_ch = 0

     # SoC dynamic
     soc_next = x_h2.soc * (1 - h2tank.η_self * Δh) - (power_ch * h2tank.η_ch + power_dch / h2tank.η_dch) * Δh / x_h2.Erated

     # State variable bounds
     overshoot = (round(soc_next;digits=3) < h2tank.α_soc_min) || (round(soc_next;digits=3) > h2tank.α_soc_max)

     overshoot ? soc_next = max(x_h2.soc * (1 - h2tank.η_self), h2tank.α_soc_min) : nothing
     overshoot ? power_ch = power_dch = 0 : nothing

     return soc_next, power_ch + power_dch
end

#                               Investment dynamic
#______________________________________________________________________________
function compute_investment_dynamics(h2tank::H2Tank, x_tank, u_tank)
     #=
         INPUT :
                 x_tank = [Erated[y], soc[end,y]]
                 u_tank[y] = h2tank control inv in kWh
         OUTPUT :
                 E_next
                 soc_next
     =#

     # Model
     if round(u_tank) > 0.
         E_next = u_tank
         soc_next = h2tank.soc[1,1,1]
     else
         E_next = x_tank.Erated
         soc_next = x_tank.soc
     end

     return E_next, soc_next
end
