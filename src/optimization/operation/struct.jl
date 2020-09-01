#=
    Controller structure for the operation
=#

# Rule based
mutable struct RuleBasedController <: AbstractController
    u::NamedTuple
    π::Function
    RuleBasedController() = new()
end

# Anticipative
mutable struct AnticipativeController <: AbstractController
    u::NamedTuple
    horizon::Int64
    model::JuMP.Model
    AnticipativeController() = new()
end

# Model predictive control
mutable struct MPCController <: AbstractController
    u::NamedTuple
    horizon::Int64
    markovchains::NamedTuple
    model::JuMP.Model
    MPCController() = new()
end
