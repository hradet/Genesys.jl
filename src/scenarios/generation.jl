#=
    Scenario generation methods
=#
abstract type AbstractScenariosGenerator end

mutable struct MarkovChain
    states
    transition_matrices
end

# Compute markov chain
function MarkovChain(data, nstate::Int64, algo::String)
    # Parameters
    nk = size(data,1)
    nh = 24
    nm = size(data[1],1)

    # Pre-allocate
    states = [zeros(nh, nk, nstate) for m in 1:nm]
    transition_matrices = [zeros(nstate, nstate, nh - 1) for m in 1:nm]

    # For each month, we compute the markov states and transition matrices
    for m in 1:nm
        # For each hour, we extract markov states using k-means from aggregated data to keep the synchronicity
        states[m], sequences = compute_markov_states([data[k][m] for k in 1:nk], nstate, algo)
        # We compute the transition matrices between consecutive hours using the sequences
        transition_matrices[m] = compute_transition_matrices_between_hours(sequences, nstate)
    end

    return MarkovChain(states, transition_matrices)
end
# Compute states using clustering algorithm
function compute_markov_states(data, nstate::Int64, algo::String)
    # Parameters
    nk = size(data,1)
    nh = size(data[1], 1)
    nd = size(data[1], 2)

    # Pre-allocate
    states = zeros(nh, nk, nstate)
    sequences = zeros(Int64, nh+1, nd)

    # Normalization
    data_n = data ./ maximum.(data)

    # Clustering
    for h in 1:nh
        if algo == "kmeans"
            # Aggregate normalized data
            data_agg = permutedims(hcat([data_n[k][h,:] for k in 1:nk]...))
            # Clustering
            clusters = kmeans(data_agg, nstate)
            # Store the denormalized states and sequences
            states[h,:,:] = clusters.centers .* maximum.(data)
            sequences[h,:] = permutedims(clusters.assignments)
        elseif algo == "kmedoids"
            # Aggregate normalized data
            data_agg = permutedims(hcat([data_n[k][h,:] for k in 1:nk]...))
            # Compute euclidean distance matrix
            dist = pairwise(Euclidean(), data_agg, dims = 2 )
            # Clustering
            clusters = kmedoids(dist, nstate)
            # Store the states and sequences
            states[h,:,:] = data_agg[:, clusters.medoids] .* maximum.(data)
            sequences[h,:] = permutedims(clusters.assignments)
        end
    end
    # The sequence of the last hour is equal to the first one
    sequences[end, :] = sequences[1, :]

    return states, sequences
end
# Compute transition matrix between consecutive hours
function compute_transition_matrices_between_hours(sequence, nstate::Int64)
    # Parameters
    nh = size(sequence, 1)

    # Pre-allocate
    transition_matrix = zeros(nstate, nstate, nh - 1)

    for h in 1:nh-1
        for (i,j) in zip(sequence[h, :], sequence[h+1, :])
            transition_matrix[i, j, h] += 1
        end

        # Convert into probabilities
        for i in 1:size(transition_matrix,1)

            sum_row = sum(transition_matrix[i, :, h])

            if sum_row > 0 # to avoid zero division...
                transition_matrix[i, :, h] = transition_matrix[i, :, h] ./ sum_row
            elseif sum_row == 0
                transition_matrix[i, i, h] = 1 # to avoid bug in Categorical...
            else
                print("Error! The sum couldn't be negative")
                return
            end
        end
    end

    return transition_matrix
end
# Clustering data by week for each month
function clustering_month_week(data, t; flag::String="")
    # Parameters
    nm = 12 # number of month in one year

    # Pre_allocate
    data_cluster= Vector{Array{Float64,2}}(undef, nm)

    for m in 1:nm
        # We cluster week days using the "not" operator
        if flag == "odd" # odd indices
            data_cluster[m] = reshape(data[(Dates.month.(t) .== m) .& .~isweekend(t)], 24, :)[:, isodd.(1:end)]
        elseif flag == "even" # even indices
            data_cluster[m] = reshape(data[(Dates.month.(t) .== m) .& .~isweekend(t)], 24, :)[:, iseven.(1:end)]
        else
            data_cluster[m] = reshape(data[(Dates.month.(t) .== m) .& .~isweekend(t)], 24, :)
        end
    end

    return data_cluster
end
# Clustering data by weekend for each month
function clustering_month_weekend(data, t; flag::String="")
    # Parameters
    nm = 12 # number of month in one year

    # Pre_allocate
    data_cluster = Vector{Array{Float64,2}}(undef, nm)

    for m in 1:nm
        # We cluster weekend days
        if flag == "odd" # odd indices
            data_cluster[m] = reshape(data[(Dates.month.(t) .== m) .& isweekend(t)], 24, :)[:, isodd.(1:end)]
        elseif flag == "even" # even indices
            data_cluster[m] = reshape(data[(Dates.month.(t) .== m) .& isweekend(t)], 24, :)[:, iseven.(1:end)]
        else
            data_cluster[m] = reshape(data[(Dates.month.(t) .== m) .& isweekend(t)], 24, :)
        end
    end

    return data_cluster
end
# Clustering data by month
function clustering_month(data, t; flag::String="")
    # Parameters
    nm = 12 # number of month in one year

    # Pre_allocate
    data_cluster = Vector{Array{Float64,2}}(undef, nm)

    for m in 1:nm
        # We cluster weekend days
        if flag == "odd" # odd indices
            data_cluster[m] = reshape(data[(Dates.month.(t) .== m)], 24, :)[:, isodd.(1:end)]
        elseif flag == "even" # even indices
            data_cluster[m] = reshape(data[(Dates.month.(t) .== m)], 24, :)[:, iseven.(1:end)]
        else
            data_cluster[m] = reshape(data[(Dates.month.(t) .== m)], 24, :)
        end
    end

    return data_cluster
end

mutable struct MarkovGenerator <: AbstractScenariosGenerator
    nstate::Int64
    algo::String
    markovchain::NamedTuple{(:wk, :wkd), Tuple{MarkovChain, MarkovChain}}

    MarkovGenerator(; nstate = 10, algo = "kmeans") = new(nstate, algo)
end

function initialize_generator!(generator::MarkovGenerator, data...)
    # Clustering by week and weekend for each month
    wk = [clustering_month_week(d.power, d.t) for d in data]
    wkd = [clustering_month_weekend(d.power, d.t) for d in data]

    # The number of states must be greater than 2 and lower than the minimum number of week and weekend days
    min_wk, min_wkd =  minimum([size(wk[k][m])[2] for k in 1:size(wk,1), m in 1:12]), minimum([size(wkd[k][m])[2] for k in 1:size(wk,1), m in 1:12])
    1 < generator.nstate < min(min_wk, min_wkd) ? nothing : generator.nstate = min(min_wk, min_wkd) - 1

    # Compute markov chains
    generator.markovchain = (wk = MarkovChain(wk, generator.nstate, generator.algo),
                             wkd = MarkovChain(wkd, generator.nstate, generator.algo))

    return generator
end

# Generate scenario from markov chains
function generate(generator::MarkovGenerator, s0, t0::DateTime, nstep::Int64; ny::Int64=1, ns::Int64=1)
    # Initialize state with the closest value
    mc = chose(generator, t0)
    Δ_s0 = [s0 .- mc.states[Dates.month(t0)][Dates.hour(t0)+1, :, :][:,state] for state in 1:generator.nstate]
    idx_0 = findmin(norm.(Δ_s0))[2]

    # Timestamp
    t = t0:Hour(1):t0+Hour(nstep-1)
    hour, month = Dates.hour.(t) .+ 1, Dates.month.(t) # hour(t[1:24]) = 0:23...

    # Pre-allocate
    ω_generated = [zeros(nstep, ny, ns) for _ in 1:size(mc.states[1],2)]
    probabilities = ones(nstep, ny, ns)

    # Initialization
    for j in 1:size(mc.states[1],2)
        ω_generated[j][1,:,:] .= s0[j]
    end

    # Generate scenarios
    for s in 1:ns, y in 1:ny, h in 1:nstep-1
        # Test to chose the appropriate markov chain
        mc = chose(generator, t[h])
        # Build categorical distribution with transition matrices - we make the assumption that the probabilities from 24h to 1h are equal whatever day transitions
        distribution = Categorical(mc.transition_matrices[month[h]][idx_0, :, hour[h]])
        # Sample state index from the distribution
        idx_1 = rand(distribution)
        # Retrieve the associated values
        for j in 1:size(mc.states[1],2)
            ω_generated[j][h+1,y,s] = mc.states[month[h]][hour[h+1], j, idx_1]
        end
        # Store probabilities
        probabilities[h+1,y,s] = probs(distribution)[idx_1]
        # Update idx
        idx_0 = idx_1
    end

    return ω_generated, prod(probabilities, dims=1)[1,:,:] / sum(prod(probabilities, dims=1)[1,:,:])
end

# Perfect foresight generator
mutable struct AnticipativeGenerator <: AbstractScenariosGenerator
    forecast
    AnticipativeGenerator() = new()
end

function initialize_generator!(generator::AnticipativeGenerator, data...)
    # Store anticipative forecast
    generator.forecast = [d for d in data]
    return generator
end

# Generate perfect forecast
function generate(generator::AnticipativeGenerator, s0, t0::DateTime, nstep::Int64; ny::Int64=1, ns::Int64=1)
    # Current index
    idx = findfirst(generator.forecast[1].t .== t0)[1]
    if length(generator.forecast[1].t[idx:end]) < nstep
        # Windows
        windows = idx : idx + length(generator.forecast[1].t[idx:end]) - 1
        # Add zeros to have a constant size
        n_zeros = nstep - length(windows)
        # Forecast
        forecast = [vcat(d.power[windows], zeros(n_zeros)) for d in generator.forecast]
    else
        # Windows
        windows = idx : idx + nstep - 1
        # Forecast
        forecast = [d.power[windows] for d in generator.forecast]
    end

    return forecast, [1.]
end
