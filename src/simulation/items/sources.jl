#=
    Sources modelling
 =#

mutable struct Source
     lifetime::Float64
     # Initial conditions
     powerMax_ini
     # Variables
     power_E::AbstractArray{Float64,3}
     powerMax::AbstractArray{Float64,2}
     # Eco
     C_pv::AbstractArray{Float64,2}
     # Inner constructor
     Source(;lifetime=25, powerMax_ini = 0.) = new(lifetime, powerMax_ini)
end

### Preallocation
function preallocate!(pv::Source, nh::Int64, ny::Int64, ns::Int64)
     pv.power_E = convert(SharedArray,zeros(nh, ny, ns))
     pv.powerMax = convert(SharedArray,zeros(ny+1, ns)) ; pv.powerMax[1,:] .= pv.powerMax_ini
     pv.C_pv = convert(SharedArray,zeros(ny, ns))
 end

 ### Operation dynamic

 ### Investment dynamic
 function compute_investment_dynamics(pv::Source, x_pv::NamedTuple, u_pv::Union{Float64, Int64})
     #=
         INPUT :
                 x_pv = [powerMax[y]]
                 u_pv[y] = pv control inv in kW
         OUTPUT :
                 pMax_next
     =#

     # Model
     if round(u_pv) > 0.
         powerMax_next = u_pv
     else
         powerMax_next = x_pv.powerMax
     end

     return powerMax_next
 end
