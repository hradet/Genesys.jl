#=
    Designer structure for the investment
=#

# Rule based
mutable struct RuleBasedDesigner <: AbstractDesigner
    u::NamedTuple
    π::Function
    RuleBasedDesigner() = new()
end

# Anticipative multi-stage
mutable struct AnticipativeMultiStageDesigner <: AbstractDesigner
    u::NamedTuple
    horizon::Int64
    model::JuMP.Model
    AnticipativeMultiStageDesigner() = new()
end

# Anticipative one-stage
mutable struct AnticipativeOneStageDesigner <: AbstractDesigner
    u::NamedTuple
    horizon::Int64
    model::JuMP.Model
    AnticipativeOneStageDesigner() = new()
end

# Anticipative one-stage with online optimization
mutable struct AnticipativeOneStageOnlineDesigner <: AbstractDesigner
    u::NamedTuple
    horizon::Int64
    model::JuMP.Model
    AnticipativeOneStageOnlineDesigner() = new()
end
