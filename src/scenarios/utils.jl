#=
    Utility functions
=#

# True if t is a weeken day
isweekend(t::Union{DateTime, Array{DateTime,1}, Array{DateTime,2}}) = (Dates.dayname.(t) .== "Saturday") .| (Dates.dayname.(t) .== "Sunday")

# Chose the right markovchain as a function of t
chose(generator::MarkovGenerator, t::DateTime) = isweekend(t) ? generator.markovchain.wkd : generator.markovchain.wk

# Generate a yearly profile from typical days clustered data
generate(data_td::Array{Float64,2}, assignments::Array{Int64,1}) = reshape(data_td[:, assignments, :], size(data_td, 1) * size(assignments, 1))
