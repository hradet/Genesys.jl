abstract type EnergyCarrier end

mutable struct Electricity <: EnergyCarrier
    in::AbstractArray{Float64,3}
    out::AbstractArray{Float64,3}

    Electricity() = new()
end

mutable struct Heat <: EnergyCarrier
    in::AbstractArray{Float64,3}
    out::AbstractArray{Float64,3}

    Heat() = new()
end

mutable struct Hydrogen <: EnergyCarrier
    in::AbstractArray{Float64,3}
    out::AbstractArray{Float64,3}

    Hydrogen() = new()
end
