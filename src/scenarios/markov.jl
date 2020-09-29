
mutable struct MarkovChain
    states
    transition_matrices
end


# Clustering functions
# Clustering data by month
function clustering_month(data, timestamp)
    #=
        Data are only clustered by month from the time array
    =#

    # Parameters
    ny = size(data,2) # number of years
    nm = 12 # number of month in one year

    # Pre_allocate
    data_cluster = Array{Array{Float64,2}}(undef, nm, ny)

    for y in 1:ny, m in 1:nm
        # We cluster data by month
        data_cluster[m,y] = reshape(data[:,y,:][Dates.month.(timestamp[:,y,:]) .== m], 24, :)
    end

    return data_cluster
end
# Clustering data by week for each month
function clustering_month_week(data, timestamp)
    #=
        Data are clustered by week days for each month from the
        time array
    =#

    # Parameters
    ny = size(data,2) # number of years
    nm = 12 # number of month in one year

    # Pre_allocate
    data_cluster= Array{Array{Float64,2}}(undef, nm, ny)

    for y in 1:ny, m in 1:nm
        # We cluster week days using the "not" operator
        data_cluster[m,y] = reshape(data[:,y,:][(Dates.month.(timestamp[:,y,:]) .== m) .& .~isweekend(timestamp[:,y,:])], 24, :)
    end

    return data_cluster
end
# Clustering data by weekend for each month
function clustering_month_weekend(data, timestamp)
    #=
        Data are clustered by week-end days for each month from the
        time array
    =#

    # Parameters
    ny = size(data,2) # number of years
    nm = 12 # number of month in one year

    # Pre_allocate
    data_cluster = Array{Array{Float64,2}}(undef, nm, ny)

    for y in 1:ny, m in 1:nm
        # We cluster weekend days
        data_cluster[m,y] = reshape(data[:,y,:][(Dates.month.(timestamp[:,y,:]) .== m) .& isweekend(timestamp[:,y,:])], 24, :)
    end

    return data_cluster
end
# Clustering states using k-means algorithm
function clustering_states(data_day, n_state)
    # Parameters
    nh = size(data_day,1)
    nd = size(data_day,2)
    # Pre-allocate
    sequence = zeros(Int64, nh, nd)
    states = zeros(nh, n_state)

    # For each hour, we extract markov states and the corresponding sequence
    for h in 1 : nh
        cluster = kmeans(convert(Array{Float64,2},reshape(data_day[h,:],1,:)), n_state)
        states[h,:], sequence[h,:] = cluster.centers, cluster.assignments
    end

    return states, sequence
end

# Transition matrix
# Compute transition matrix between consecutive hours
function compute_transition_between_hours(sequence, n_state::Int64)
    # Parameters
    nh = size(sequence, 1)

    # Pre-allocate
    transition_matrix = zeros(n_state,n_state,nh-1)

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
# Compute transition matrix from a sequence
function compute_transition_from_sequence(sequence, n_state::Int64)

    transition_matrix = zeros(n_state,n_state) # n of state size square matrix

    # Compute number of transition
    for (i,j) in zip(sequence[1:end-1],sequence[2:end])
        transition_matrix[i,j] += 1
    end

    # Convert into probabilities
    for i in 1:nMarkovstate
        sum_row = sum(transition_matrix[i,:])

        if sum_row > 0 # to avoid zero division...
            transition_matrix[i,:] = transition_matrix[i,:] ./ sum_row
        elseif sum_row == 0
            transition_matrix[i, i] = 1 # to avoid bug in Categorical...
        else
            print("Error! The sum couldn't be negative")
            return
        end
    end

    return transition_matrix
end

# Compute MarkovChain
# Compute markov-chain from clustered data
function compute_markovchain(data, n_state::Int64)
    # Parameters
    nm = size(data,1)
    ny = size(data,2)

    # Pre-allocate
    states = Array{Array{Float64,2}}(undef, nm, ny)
    transition_matrices = Array{Array{Float64,3}}(undef, nm, ny)

    for y in 1:ny
        # Compute transition matrices and states for each month
        for m in 1:nm
            # # For each hour, we extract markov states using k-means and the sequence
            # needed to compute transition matrices between consecutive hours. The number
            # of markostate must be < number of days...
            states[m,y], sequence = clustering_states(data[m,y], min(n_state, size(data[m,y],2)-1))

            # Compute transition matrix between consecutive hours for each month
            transition_matrices[m,y] = compute_transition_between_hours(sequence, min(n_state, size(data[m,y],2)-1))
        end
    end

    return MarkovChain(states, transition_matrices)
end

# Compute scenario from markov_chain
# Compute a markov scenario
function compute_scenario(mc::MarkovChain, state_init::Union{Int64,Float64}, time_init::DateTime, year::Int64, nstep::Int64)
    # Initialize timestamp
    timestamp = time_init:Hour(1):time_init+Hour(nstep)

    # Initialize scenario with current value
    scenario = zeros(length(timestamp))
    scenario[1] = state_init

    # Initialize state with the closest value
    idx_states = zeros(Int64, length(timestamp))
    idx_states[1] = findmin(abs.(state_init .- mc.states[Dates.month(time_init), year][Dates.hour(time_init)+1, :]))[2] # hour starts at 0...

    # Simulate
    for k in 2:length(timestamp)
        # Retrieve the current hour and month
        hour, month = Dates.hour(timestamp[k])+1, Dates.month(timestamp[k])
        # Build categorical distribution with transition matrices
        distributions = [Categorical(mc.transition_matrices[month, year][s, :, max(1,hour-1)]) for s in 1:size(mc.transition_matrices[month,year],1)]
        # Compute state index from the probability matrix
        idx_states[k] = rand(distributions[idx_states[k-1]])
        # Retrieve the associated values
        scenario[k] = mc.states[month, year][hour, idx_states[k]]
    end

    return scenario
end
# Compute a markov scenario with multiple markov chain
function compute_scenario(mc_wk::MarkovChain, mc_wkd::MarkovChain, state_init::Union{Int64,Float64}, time_init::DateTime, year::Int64, nstep::Int64)
    # Initialize timestamp
    timestamp = time_init:Hour(1):time_init+Hour(nstep)

    # Initialize scenario with current value
    scenario = zeros(length(timestamp))
    scenario[1] = state_init

    # Initialize state with the closest value
    idx_states = zeros(Int64, length(timestamp))
    isweekend(timestamp[1]) ? mc = mc_wkd : mc = mc_wk # test to chose the appropriate markov chain
    idx_states[1] = findmin(abs.(state_init .- mc.states[Dates.month(time_init), year][Dates.hour(time_init)+1, :]))[2]

    # Simulate
    for k in 2:length(timestamp)
        # Retrieve the current hour and month
        hour, month = Dates.hour(timestamp[k])+1, Dates.month(timestamp[k])
        # Test to chose the appropriate markov chain
        isweekend(timestamp[k]) ? mc = mc_wkd : mc = mc_wk
        # Build categorical distribution with transition matrices
        distributions = [Categorical(mc.transition_matrices[month, year][s, :, max(1,hour-1)]) for s in 1:size(mc.transition_matrices[month,year],1)]
        # Compute state index from the probability matrix
        idx_states[k] = rand(distributions[idx_states[k-1]])
        # Retrieve the associated values
        scenario[k] = mc.states[month, year][hour, idx_states[k]]
    end

    return scenario
end
# Compute scenarios from simple sequence
function compute_scenario(mc::MarkovChain, state_init, horizon::Int64)
    # init = NaN => initial state randomly chosen
    @assert size(mc.transition_matrices)[1] == size(mc.transition_matrices)[2] # square required

    # create vector of discrete RVs for each row
    distributions = [Categorical(mc.transition_matrices[i, :]) for i in 1:size(mc.transition_matrices)[1]]

    # setup the simulation
    scenario = zeros(Int64, horizon)
    scenario[1] = state_init # set the initial state

    for k in 2:horizon
        scenario[k] = rand(distributions[scenario[k-1]]) # draw new value from last state's transition distribution
    end

    return scenario
end

# Compute markovchains from scenarios
function compute_markovchains(ω_optim::Scenarios)
    # TODO : ajouter struct avec mc, nscenarios, nstate
    # Clustering data
    pv = clustering_month(ω_optim.values.pv_E, ω_optim.timestamp)
    ld_E_wk = clustering_month_week(ω_optim.values.ld_E, ω_optim.timestamp)
    ld_E_wkd = clustering_month_weekend(ω_optim.values.ld_E, ω_optim.timestamp)
    ld_H_wk = clustering_month_week(ω_optim.values.ld_H, ω_optim.timestamp)
    ld_H_wkd = clustering_month_weekend(ω_optim.values.ld_H, ω_optim.timestamp)
    # Compute markov chain
    n_state = minimum(vcat(size.(pv[:,1],2), size.(ld_E_wk[:,1],2), size.(ld_E_wkd[:,1],2)))-1
    mc_pv = compute_markovchain(pv, n_state)
    mc_ld_E_wk = compute_markovchain(ld_E_wk, n_state)
    mc_ld_E_wkd = compute_markovchain(ld_E_wkd, n_state)
    mc_ld_H_wk = compute_markovchain(ld_H_wk, n_state)
    mc_ld_H_wkd = compute_markovchain(ld_H_wkd, n_state)
    # Store in a Namedtuple
    markovchains = (
    pv_E = mc_pv,
    ld_E = (wk = mc_ld_E_wk, wkd = mc_ld_E_wkd),
    ld_H = (wk = mc_ld_H_wk, wkd = mc_ld_H_wkd),
    )

    return markovchains
end
