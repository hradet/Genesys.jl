#=
    Li-ion battery modelling
 =#

 #                                  Structure
 #______________________________________________________________________________
 struct Liion
     # Paramètres
     α_p_ch::Float64
     α_p_dch::Float64
     η_ch::Float64
     η_dch::Float64
     η_self::Float64
     α_soc_min::Float64
     α_soc_max::Float64
     lifetime::Float64
     nCycle::Float64
     dod::Float64
     # Variables
     Erated::AbstractArray{Float64,2}
     power_E::AbstractArray{Float64,3}
     soc::AbstractArray{Float64,3}
     soh::AbstractArray{Float64,3}
     # Eco
     C_liion::AbstractArray{Float64,2}
 end
 # Constructor
 function Liion(outputGUI::NamedTuple, nh::Int64, ny::Int64, ns::Int64)
     # Paramètres
     α_p_ch = outputGUI.α_p_ch
     α_p_dch = outputGUI.α_p_dch
     η_ch = outputGUI.η_ch
     η_dch = outputGUI.η_dch
     η_self = outputGUI.η_self
     α_soc_min = outputGUI.α_soc_min
     α_soc_max = outputGUI.α_soc_max
     lifetime = outputGUI.lifetime
     nCycle = outputGUI.nCycle
     dod = outputGUI.dod
     # Variables
     Erated = convert(SharedArray,zeros(ny+1, ns))
     power_E = convert(SharedArray,zeros(nh, ny, ns))
     soc = convert(SharedArray,zeros(nh+1, ny+1, ns))
     soh = convert(SharedArray,zeros(nh+1, ny+1, ns))
     # Eco
     C_liion = convert(SharedArray,zeros(ny, ns))
     return Liion(α_p_ch,α_p_dch,η_ch,η_dch,η_self,α_soc_min,α_soc_max,lifetime,nCycle,dod,Erated,power_E,soc,soh,C_liion)
 end

 #                               Operation dynamic
 #______________________________________________________________________________
 function compute_operation_dynamics(liion::Liion, x_liion::NamedTuple, u_liion::Float64, Δh::Int64)
     #=
     INPUT :
             x_liion = (Erated[y], soc[h,y], soh[h,y]) tuple
             u_liion[h,y] = control power in kW
     OUTPUT :
             soc_next
             soh_next
             power = the real battery power in kW
     =#

     # Power constraint and correction
     0 <= u_liion <= liion.α_p_dch * x_liion.Erated ? power_dch = u_liion : power_dch = 0
     0 <= -u_liion <= liion.α_p_ch * x_liion.Erated ? power_ch = u_liion : power_ch = 0

     # SoC dynamic
     soc_next = x_liion.soc * (1 - liion.η_self * Δh) - (power_ch * liion.η_ch + power_dch / liion.η_dch) * Δh / x_liion.Erated

     # SoH dynamic
     soh_next = x_liion.soh - (power_dch - power_ch) * Δh / (2 * liion.nCycle * liion.dod * x_liion.Erated)

     # SoC and SoH constraints and corrections
     overshoot = (round(soc_next;digits=3) < liion.α_soc_min) || (round(soc_next;digits=3) > liion.α_soc_max) || (soh_next < 0)

     overshoot ? soc_next = max(x_liion.soc * (1 - liion.η_self), liion.α_soc_min) : nothing
     overshoot ? soh_next = x_liion.soh : nothing
     overshoot ? power_ch = power_dch = 0 : nothing

     return soc_next, soh_next, power_ch + power_dch
 end

 #                               Investment dynamic
 #______________________________________________________________________________
 function compute_investment_dynamics(liion::Liion, x_liion::NamedTuple, u_liion::Float64)
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
     if round(u_liion) > 0.
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
