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
