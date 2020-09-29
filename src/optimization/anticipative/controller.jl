# Anticipative controller

mutable struct Anticipative <: AbstractController
    u::NamedTuple
    model::JuMP.Model
    Anticipative() = new()
end

### Online
function compute_operation_decisions!(h::Int64, y::Int64, s::Int64, des::DistributedEnergySystem, controller::Anticipative)
    return controller
end
