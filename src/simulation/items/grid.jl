#=
    Elec grid modelling
 =#

mutable struct Grid
     # Parameters
     powerMax
     # Variables
     power_E::AbstractArray{Float64,3}
     C_grid_in::AbstractArray{Float64,3}
     C_grid_out::AbstractArray{Float64,3}
     # Inner constructor
     Grid(; powerMax = 36.) = new(powerMax)
end

### Preallocation
function preallocate!(grid::Grid, nh::Int64, ny::Int64, ns::Int64)
     grid.power_E = convert(SharedArray,zeros(nh, ny, ns))
     grid.C_grid_in = convert(SharedArray,zeros(nh, ny, ns))
     grid.C_grid_out = convert(SharedArray,zeros(nh, ny, ns))
end

### Operation dynamic

### Investment dynamic
