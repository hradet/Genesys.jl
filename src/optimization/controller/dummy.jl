#=
    Dummy controller
=#

mutable struct DummyController <: AbstractController
    u::NamedTuple
    DummyController() = new()
end

### Offline
function initialize_controller!(des::DistributedEnergySystem, controller::DummyController, ω::AbstractScenarios)
    # Preallocation
    preallocate!(controller, des.parameters.nh, des.parameters.ny, des.parameters.ns)
    return controller
end

### Online
function compute_operation_decisions!(h::Int64, y::Int64, s::Int64, des::DistributedEnergySystem, controller::DummyController)
    return controller
end
