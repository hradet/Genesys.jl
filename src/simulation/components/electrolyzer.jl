#=
    Electrolyzer modelling
 =#

#                                  Structure
#______________________________________________________________________________
struct Electrolyzer
     # Paramètres
     α_p::Float64
     η_E_H2::Float64
     η_E_H::Float64
     lifetime::Float64
     nHoursMax::Float64
     # Variables
     powerMax::AbstractArray{Float64,2}
     power_E::AbstractArray{Float64,3}
     power_H::AbstractArray{Float64,3}
     power_H2::AbstractArray{Float64,3}
     soh::AbstractArray{Float64,3}
     # Eco
     C_elyz::AbstractArray{Float64,2}
end
# Constructor
function Electrolyzer(outputGUI::NamedTuple, nh::Int64, ny::Int64, ns::Int64)
     # Paramètres
     α_p = outputGUI.α_p
     η_E_H2 = outputGUI.η_E_H2
     η_E_H = outputGUI.η_E_H
     lifetime = outputGUI.lifetime
     nHoursMax = outputGUI.nHoursMax
     # Variables
     powerMax = convert(SharedArray,zeros(ny+1, ns))
     power_E = convert(SharedArray,zeros(nh, ny, ns))
     power_H = convert(SharedArray,zeros(nh, ny, ns))
     power_H2 = convert(SharedArray,zeros(nh, ny, ns))
     soh = convert(SharedArray,zeros(nh+1, ny+1, ns))
     # Eco
     C_elyz = convert(SharedArray,zeros(ny, ns))
     return Electrolyzer(α_p,η_E_H2,η_E_H,lifetime,nHoursMax,powerMax,power_E,power_H,power_H2,soh,C_elyz)
end

#                               Operation dynamic
#______________________________________________________________________________
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
    elyz.α_p * x_elyz.powerMax <= -u_elyz <= x_elyz.powerMax ? power_E = u_elyz : power_E = 0

    # Power computations
    power_H2 = - power_E * elyz.η_E_H2
    power_H = - power_E * elyz.η_E_H

    # SoH computation
    soh_next = x_elyz.soh - (power_E < 0.) * Δh / elyz.nHoursMax

    # SoH constraint and correction
    soh_next < 0. ? power_E = power_H = power_H2 = 0 : nothing
    soh_next < 0. ? soh_next = x_elyz.soh : nothing

    return power_E, power_H, power_H2, soh_next
end

#                               Investment dynamic
#______________________________________________________________________________
function compute_investment_dynamics(elyz::Electrolyzer, x_elyz::NamedTuple, u_elyz::Float64)
    #=
        INPUT :
                x_elyz = [powerMax[y], soh[end,y]]
                u_elyz[y] = elyz control inv in kW
        OUTPUT :
                pMax_next
                soh_next
    =#

    # Model
    if round(u_elyz) > 0.
        powerMax_next = u_elyz
        soh_next = 1.
    else
        powerMax_next = x_elyz.powerMax
        soh_next = x_elyz.soh
    end

    return powerMax_next, soh_next
end
