mutable struct MarkovChain
    states
    transition_matrices
end

mutable struct MarkovChainOptions
    algorithm::String
    selection::String

    MarkovChainOptions(; algorithm = "kmeans", selection = "") = new(algorithm, selection)
end

# True if t is a weeken day
isweekend(t::Union{DateTime, Array{DateTime,1}, Array{DateTime,2}}) = (Dates.dayname.(t) .== "Saturday") .| (Dates.dayname.(t) .== "Sunday")

# Compute MarkovChain
# Compute markov-chain from aggregated data
function compute_markovchains(pv, ld_E, ld_H; n_state::Int64=0, options::MarkovChainOptions = MarkovChainOptions())
    # Clustering by week and weekend for each month
    pv_wk, pv_wkd = clustering_month_week(pv.power, pv.t, options.selection), clustering_month_weekend(pv.power, pv.t, options.selection)
    ld_E_wk, ld_E_wkd = clustering_month_week(ld_E.power, ld_E.t, options.selection), clustering_month_weekend(ld_E.power, ld_E.t, options.selection)
    ld_H_wk, ld_H_wkd = clustering_month_week(ld_H.power, ld_H.t, options.selection), clustering_month_weekend(ld_H.power, ld_H.t, options.selection)

    # The number of states must be greater than 2 and lower than the minimum of week and weekend days
    1 < n_state < min(minimum(size.(pv_wk,2)), minimum(size.(pv_wkd,2))) ? nothing : n_state = min(minimum(size.(pv_wk,2)), minimum(size.(pv_wkd,2))) - 1

    # Compute markov chains
    mc_wk = compute_markovchain(pv_wk, ld_E_wk, ld_H_wk, n_state, options)
    mc_wkd = compute_markovchain(pv_wkd, ld_E_wkd, ld_H_wkd, n_state, options)

    return (wkd = mc_wkd, wk = mc_wk)
end

# Compute markov-chain from aggregated data
function compute_markovchain(pv, ld_E, ld_H, n_state::Int64, options::MarkovChainOptions)
    # Parameters
    nh = 24
    nm = size(pv,1)

    # Pre-allocate
    states = [zeros(nh, 3, n_state) for m in 1:nm]
    transition_matrices = [zeros(n_state, n_state, nh - 1) for m in 1:nm]

    # For each month, we compute the markov states and transition matrices
    for m in 1:nm
        # For each hour, we extract markov states using k-means from aggregated data to keep the synchronicity
        states[m], sequences = compute_markov_states(pv[m], ld_E[m], ld_H[m], n_state, options)
        # We compute the transition matrices between consecutive hours using the sequences
        transition_matrices[m] = compute_transition_matrices_between_hours(sequences, n_state)
    end

    return MarkovChain(states, transition_matrices)
end

# Clustering states using k-means algorithm
function compute_markov_states(pv, ld_E, ld_H, n_state::Int64, options::MarkovChainOptions)
    # Parameters
    nh = size(pv, 1)
    nd = size(pv, 2)

    # Pre-allocate
    states = zeros(nh, 3, n_state)
    sequences = zeros(Int64, nh, nd)

    # Normalization
    pv_n, ld_E_n, ld_H_n = pv ./ maximum(pv), ld_E ./ maximum(ld_E), ld_H ./ maximum(ld_H)

    # Clustering
    for h in 1:nh
        if options.algorithm == "kmeans"
            # Aggregate normalized data
            data_agg = permutedims(hcat(pv_n[h,:], ld_E_n[h,:], ld_H_n[h,:]))
            # Clustering
            clusters = kmeans(data_agg, n_state)
            # Store the states and sequences
            states[h,:,:] = clusters.centers
            sequences[h,:] = permutedims(clusters.assignments)
        elseif options.algorithm == "kmeans_minmax"
            # Aggregate normalized data
            data_agg = permutedims(hcat(pv_n[h,:], ld_E_n[h,:], ld_H_n[h,:]))
            # Clustering
            clusters = kmeans(data_agg, n_state)
            # Replace min max centers by the min max real values
            clusters.centers[argmin(clusters.centers, dims=2)] = minimum(data_agg, dims=2)
            clusters.centers[argmax(clusters.centers, dims=2)] = maximum(data_agg, dims=2)
            # Store the states and sequences
            states[h,:,:] = clusters.centers
            sequences[h,:] = permutedims(clusters.assignments)
        elseif options.algorithm == "kmedoids"
            # Aggregate normalized data
            data_agg = permutedims(hcat(pv_n[h,:], ld_E_n[h,:], ld_H_n[h,:]))
            # Compute euclidean distance matrix
            dist = pairwise(Euclidean(), data_agg, dims = 2 )
            # Clustering
            clusters = kmedoids(dist, n_state)
            # Store the states and sequences
            states[h,:,:] = data_agg[:, clusters.medoids]
            sequences[h,:] = permutedims(clusters.assignments)
        end
    end

    # Denormalization
    states[:,1,:] .*= maximum(pv)
    states[:,2,:] .*= maximum(ld_E)
    states[:,3,:] .*= maximum(ld_H)

    return states, sequences
end

# Transition matrix
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
        if flag == "optim" # odd indices
            data_cluster[m] = reshape(data[(Dates.month.(t) .== m) .& .~isweekend(t)], 24, :)[:, isodd.(1:end)]
        elseif flag == "simu" # even indices
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
        if flag == "optim" # odd indices
            data_cluster[m] = reshape(data[(Dates.month.(t) .== m) .& isweekend(t)], 24, :)[:, isodd.(1:end)]
        elseif flag == "simu" # even indices
            data_cluster[m] = reshape(data[(Dates.month.(t) .== m) .& isweekend(t)], 24, :)[:, iseven.(1:end)]
        else
            data_cluster[m] = reshape(data[(Dates.month.(t) .== m) .& isweekend(t)], 24, :)
        end
    end

    return data_cluster
end

# Compute scenario from markov chains
function compute_scenarios(markovchains, s0, t0::DateTime, nstep::Int64; ny::Int64=1, ns::Int64=1)
    # Parameters
    n_state = size(markovchains.wk.transition_matrices[1],1)
    t = t0:Hour(1):t0+Hour(nstep-1)

    # Pre-allocate
    ω_pv = (t = repeat(t, 1, ny, ns),
            power = zeros(nstep, ny, ns))
    ω_ld_E = (t = repeat(t, 1, ny, ns),
             power = zeros(nstep, ny, ns))
    ω_ld_H = (t = repeat(t, 1, ny, ns),
             power = zeros(nstep, ny, ns))
    idx_s = zeros(Int64, nstep)

    # Initialization
    ω_pv.power[1,:,:] .= s0[1]
    ω_ld_E.power[1,:,:] .= s0[2]
    ω_ld_H.power[1,:,:] .= s0[3]

    # Initialize state with the closest value
    isweekend(t0) ? mc = markovchains.wkd : mc = markovchains.wk
    Δ_s0 = [s0 .- mc.states[Dates.month(t0)][Dates.hour(t0)+1, :, :][:,state] for state in 1:n_state]
    idx_s[1] = findmin(norm.(Δ_s0))[2]

    # Simulate
    for s in 1:ns, y in 1:ny, k in 2:nstep
        # Retrieve the current hour and month - hour(t[1]) = 0...
        h, m = Dates.hour(t[k]) + 1, Dates.month(t[k])
        # Test to chose the appropriate markov chain
        isweekend(t[k]) ? mc = markovchains.wkd : mc = markovchains.wk
        # Build categorical distribution with transition matrices
        distributions = [Categorical(mc.transition_matrices[m][state, :, max(1, h - 1)]) for state in 1:n_state]
        # Compute state index from the probability matrix
        idx_s[k] = rand(distributions[idx_s[k-1]])
        # Retrieve the associated values
        ω_pv.power[k,y,s] = mc.states[m][h, 1, idx_s[k]]
        ω_ld_E.power[k,y,s] = mc.states[m][h, 2, idx_s[k]]
        ω_ld_H.power[k,y,s] = mc.states[m][h, 3, idx_s[k]]
    end

    return ω_pv, ω_ld_E, ω_ld_H
end
