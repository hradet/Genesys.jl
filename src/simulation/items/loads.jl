#=
    Loads modelling
 =#

mutable struct Load
     # Variables
     power::AbstractArray{Float64,3}
     # Inner constructor
     Load() = new()
end

### Preallocation
function preallocate!(ld::Load, nh::Int64, ny::Int64, ns::Int64)
     ld.power = convert(SharedArray,zeros(nh, ny, ns))
end

 ### Operation dynamic

 ### Investment dynamic
