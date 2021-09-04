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
    pv.carrier.out = convert(SharedArray,zeros(nh, ny, ns))
    pv.powerMax = convert(SharedArray,zeros(ny+1, ns)) ; pv.powerMax[1,:] .= pv.powerMax_ini
    pv.timestamp = Array{DateTime}(undef,(nh, ny, ns))
    pv.cost = convert(SharedArray,zeros(ny, ns))

    return pv
 end

 ### Investment dynamic
 function compute_investment_dynamics!(y::Int64, s::Int64, pv::Solar, decision::Union{Float64, Int64})
     if decision > 1e-2
         pv.powerMax[y+1,s] = decision
     else
         pv.powerMax[y+1,s] = pv.powerMax[y,s]
     end
 end
