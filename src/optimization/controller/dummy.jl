#=
    Dummy controller
=#

mutable struct DummyController <: AbstractController
    u::NamedTuple
    DummyController() = new()
end

### Offline
function initialize_controller!(des::DES, ω::Scenarios)
    return nothing
end

### Online 
function compute_operation_decisions!(h::Int64, y::Int64, s::Int64, des::DES)
    return nothing
end
