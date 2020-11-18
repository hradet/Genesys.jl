mutable struct MarkovChain
    states
    transition_matrices
end

mutable struct MarkovChainOptions
    days_selection::String
    clustering_days::String
    clustering_state_algorithm::String

    MarkovChainOptions(; days_selection = "", clustering_days = "wk_wkd", clustering_state_algorithm = "kmeans") = new(days_selection, clustering_days, clustering_state_algorithm)
end

# True if t is a weeken day
isweekend(t::Union{DateTime, Array{DateTime,1}, Array{DateTime,2}}) = (Dates.dayname.(t) .== "Saturday") .| (Dates.dayname.(t) .== "Sunday")

# Compute week and wkd markov chains from aggregated input data
function compute_markovchains(data...; n_state::Int64=0, options::MarkovChainOptions = MarkovChainOptions())
    if options.clustering_days == "wk_wkd"
        # Clustering by week and weekend for each month
        wk = [clustering_month_week(d.power, d.t, options.clustering_days) for d in data]
        wkd = [clustering_month_weekend(d.power, d.t, options.clustering_days) for d in data]

        # The number of states must be greater than 2 and lower than the minimum number of week and weekend days
        min_wk, min_wkd =  minimum([size(wk[k][m])[2] for k in 1:size(wk,1), m in 1:12]), minimum([size(wkd[k][m])[2] for k in 1:size(wk,1), m in 1:12])
        1 < n_state < min(min_wk, min_wkd) ? nothing : n_state = min(min_wk, min_wkd) - 1

        # Compute markov chains
        mc_wk = compute_markovchain(wk, n_state, options)
        mc_wkd = compute_markovchain(wkd, n_state, options)

        return mc_wkd, mc_wk
    elseif options.clustering_days == "month"
        # Clustering by month
        month = [clustering_month(d.power, d.t, options.clustering_days) for d in data]

        # The number of states must be greater than 2 and lower than the minimum number of month days
        1 < n_state < minimum([size(month[k][m])[2] for k in 1:size(month,1), m in 1:12]) ? nothing : n_state = minimum([size(month[k][m])[2] for k in 1:size(month,1), m in 1:12]) - 1

        # Compute markov chains
        mc_month = compute_markovchain(month, n_state, options)

        return mc_month
    else
        prinln("Clustering day type unknown...")
    end
end
# Compute markov chain
function compute_markovchain(data, n_state::Int64, options::MarkovChainOptions)
    # Parameters
    nk = size(data,1)
    nh = 24
    nm = size(data[1],1)

    # Pre-allocate
    states = [zeros(nh, nk, n_state) for m in 1:nm]
    transition_matrices = [zeros(n_state, n_state, nh - 1) for m in 1:nm]

    # For each month, we compute the markov states and transition matrices
    for m in 1:nm
        # For each hour, we extract markov states using k-means from aggregated data to keep the synchronicity
        states[m], sequences = compute_markov_states([data[k][m] for k in 1:nk], n_state, options)
        # We compute the transition matrices between consecutive hours using the sequences
        transition_matrices[m] = compute_transition_matrices_between_hours(sequences, n_state)
    end

    return MarkovChain(states, transition_matrices)
end
# Compute states using clustering algorithm
function compute_markov_states(data, n_state::Int64, options::MarkovChainOptions)
    # Parameters
    nk = size(data,1)
    nh = size(data[1], 1)
    nd = size(data[1], 2)

    # Pre-allocate
    states = zeros(nh, nk, n_state)
    sequences = zeros(Int64, nh, nd)

    # Normalization
    data_n = data ./ maximum.(data)

    # Clustering
    for h in 1:nh
        if options.clustering_state_algorithm == "kmeans"
            # Aggregate normalized data
            data_agg = permutedims(hcat([data_n[k][h,:] for k in 1:nk]...))
            # Clustering
            clusters = kmeans(data_agg, n_state)
            # Store the denormalized states and sequences
            states[h,:,:] = clusters.centers .* maximum.(data)
            sequences[h,:] = permutedims(clusters.assignments)
        elseif options.clustering_state_algorithm == "kmedoids"
            # Aggregate normalized data
            data_agg = permutedims(hcat([data_n[k][h,:] for k in 1:nk]...))
            # Compute euclidean distance matrix
            dist = pairwise(Euclidean(), data_agg, dims = 2 )
            # Clustering
            clusters = kmedoids(dist, n_state)
            # Store the states and sequences
            states[h,:,:] = data_agg[:, clusters.medoids] .* maximum.(data)
            sequences[h,:] = permutedims(clusters.assignments)
        end
    end

    return states, sequences
end
# Compute transition matrix between consecutive hours
function compute_transition_matrices_between_hours(sequence, n_state::Int64)
    # Parameters
    nh = size(sequence, 1)

    # Pre-allocate
    transition_matrix = zeros(n_state, n_state, nh - 1)

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
function clustering_month_week(data, t, flag::String)
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
function clustering_month_weekend(data, t, flag::String)
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
function clustering_month(data, t, flag::String)
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
# Compute scenario from markov chains
function compute_scenarios(mc_wk, mc_wkd, s0, t0::DateTime, nstep::Int64; ny::Int64=1, ns::Int64=1)
    # Parameters
    n = size(mc_wk.states[1],2)
    n_state = size(mc_wk.transition_matrices[1],1)
    t = t0:Hour(1):t0+Hour(nstep-1)

    # Pre-allocate
    scenarios = [(t = repeat(t, 1, ny, ns), power = zeros(nstep, ny, ns)) for _ in 1:n]

    # Initialization
    for j in 1:n
        scenarios[j].power[1,:,:] .= s0[j]
    end

    # Initialize state with the closest value
    isweekend(t0) ? mc = mc_wkd : mc = mc_wk
    Δ_s0 = [s0 .- mc.states[Dates.month(t0)][Dates.hour(t0)+1, :, :][:,state] for state in 1:n_state]
    idx_prev = findmin(norm.(Δ_s0))[2]

    # Simulate
    for s in 1:ns, y in 1:ny, k in 2:nstep
        # Retrieve the current hour and month - hour(t[1]) = 0...
        h, m = Dates.hour(t[k]) + 1, Dates.month(t[k])
        # Test to chose the appropriate markov chain
        isweekend(t[k]) ? mc = mc_wkd : mc = mc_wk
        # Build categorical distribution with transition matrices
        distributions = [Categorical(mc.transition_matrices[m][state, :, max(1, h - 1)]) for state in 1:n_state]
        # Compute state index from the probability matrix
        idx = rand(distributions[idx_prev])
        # Retrieve the associated values
        for j in 1:n
            scenarios[j].power[k,y,s] = mc.states[m][h, j, idx]
        end
        # Update idx
        idx_prev = idx
    end

    return scenarios
end
# Compute scenario from a single markov chain
function compute_scenarios(mc, s0, t0::DateTime, nstep::Int64; ny::Int64=1, ns::Int64=1)
    # Parameters
    n = size(mc.states[1],2)
    n_state = size(mc.transition_matrices[1],1)
    t = t0:Hour(1):t0+Hour(nstep-1)

    # Pre-allocate
    scenarios = [(t = repeat(t, 1, ny, ns), power = zeros(nstep, ny, ns)) for _ in 1:n]
    idx_s = zeros(Int64, nstep)

    # Initialization
    for j in 1:n
        scenarios[j].power[1,:,:] .= s0[j]
    end

    # Initialize state with the closest value
    Δ_s0 = [s0 .- mc.states[Dates.month(t0)][Dates.hour(t0)+1, :, :][:,state] for state in 1:n_state]
    idx_s[1] = findmin(norm.(Δ_s0))[2]

    # Simulate
    for s in 1:ns, y in 1:ny, k in 2:nstep
        # Retrieve the current hour and month - hour(t[1]) = 0...
        h, m = Dates.hour(t[k]) + 1, Dates.month(t[k])
        # Build categorical distribution with transition matrices
        distributions = [Categorical(mc.transition_matrices[m][state, :, max(1, h - 1)]) for state in 1:n_state]
        # Compute state index from the probability matrix
        idx_s[k] = rand(distributions[idx_s[k-1]])
        # Retrieve the associated values
        for j in 1:n
            scenarios[j].power[k,y,s] = mc.states[m][h, j, idx_s[k]]
        end
    end

    return scenarios
end
