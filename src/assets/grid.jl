#=
    Elec grid modelling
 =#

mutable struct Grid
     # Parameters
     powerMax
     # Variables
     power_E::AbstractArray{Float64,3}
     cost_in::AbstractArray{Float64,3}
     cost_out::AbstractArray{Float64,3}
     # Inner constructor
     Grid(; powerMax = 36.) = new(powerMax)
end

### Preallocation
function preallocate!(grid::Grid, nh::Int64, ny::Int64, ns::Int64)
     grid.power_E = convert(SharedArray,zeros(nh, ny, ns))
     grid.cost_in = convert(SharedArray,zeros(nh, ny, ns))
     grid.cost_out = convert(SharedArray,zeros(nh, ny, ns))
end

### Operation dynamic

### Investment dynamic
