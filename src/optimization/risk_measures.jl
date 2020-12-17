#=
    Risk measures for optimization
=#

abstract type AbstractRiskMeasure end

# Expectation
struct Expectation <: AbstractRiskMeasure end

# Worst case
struct WorstCase <: AbstractRiskMeasure end

# Conditional value at risk
struct CVaR <: AbstractRiskMeasure
    β::Float64

    function CVaR(β::Float64)
        @assert 0. <= β <= 1. "β must be in [0,1]"
        return new(β)
    end
end

function conditional_value_at_risk(support::Array{Float64,1}, probabilities::Array{Float64,1}, α::Float64)
    # Value at risk
    var = value_at_risk(support, probabilities, α)
    # Tail values
    tail = support .< var
    return (sum(probabilities[tail] .* support[tail]) - (sum(probabilities[tail]) - α) * var) / α
end

function value_at_risk(support::Array{Float64,1}, probabilities::Array{Float64,1}, α::Float64)
    i = findfirst(cumsum(probabilities[sortperm(support)]) .>= α)
    return sort(support)[i]
end
