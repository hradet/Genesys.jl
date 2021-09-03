#=
    Elec grid modelling
 =#

mutable struct Grid <: AbstractGrid
    # Parameters
    powerMax::Float64
    # Variables
    carrier::EnergyCarrier
    cost_in::AbstractArray{Float64,3}
    cost_out::AbstractArray{Float64,3}
    # Inner constructor
    Grid(; powerMax = 36., carrier = Electricity()) = new(powerMax, carrier)
end

### Preallocation
function preallocate!(grid::Grid, nh::Int64, ny::Int64, ns::Int64)
     grid.carrier.in = convert(SharedArray,zeros(nh, ny, ns))
     grid.carrier.out = convert(SharedArray,zeros(nh, ny, ns))
     grid.cost_in = convert(SharedArray,zeros(nh, ny, ns))
     grid.cost_out = convert(SharedArray,zeros(nh, ny, ns))
     return grid
end
