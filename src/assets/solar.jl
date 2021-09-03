#=
    Sources modelling
 =#
mutable struct Solar <: AbstractGeneration
     lifetime::Int64
     # Initial conditions
     powerMax_ini::Float64
     soh_ini::Float64
     # Variables
     carrier::Electricity
     powerMax::AbstractArray{Float64,2}
     timestamp::Array{DateTime,3}
     # Eco
     cost::AbstractArray{Float64,2}
     # Inner constructor
     Solar(;lifetime=25, powerMax_ini = 0., soh_ini = 1.) = new(lifetime, powerMax_ini, soh_ini)
end

### Preallocation
function preallocate!(pv::Solar, nh::Int64, ny::Int64, ns::Int64)
    pv.carrier = Electricity()
    pv.carrier.in = convert(SharedArray,zeros(nh, ny, ns))
    pv.powerMax = convert(SharedArray,zeros(ny+1, ns)) ; pv.powerMax[1,:] .= pv.powerMax_ini
    pv.timestamp = Array{DateTime}(undef,(nh, ny, ns))
    pv.cost = convert(SharedArray,zeros(ny, ns))

    return pv
 end

 ### Operation dynamic

 ### Investment dynamic
 function compute_investment_dynamics(pv::Solar, x_pv::NamedTuple{(:powerMax,), Tuple{Float64}}, u_pv::Union{Float64, Int64})
     #=
         INPUT :
                 x_pv = [powerMax[y]]
                 u_pv[y] = pv control inv in kW
         OUTPUT :
                 pMax_next
     =#

     # Model
     if u_pv > 1e-2
         powerMax_next = u_pv
     else
         powerMax_next = x_pv.powerMax
     end

     return powerMax_next
 end
