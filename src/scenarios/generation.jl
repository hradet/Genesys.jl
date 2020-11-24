#=
    Scenario generation methods
=#


mutable struct MarkovGeneration


end

function generate(generator::MarkovGeneration)

    return ω_generated, proba
end

# Generate scenario from markov chains
function generate(mc_wk::MarkovChain, mc_wkd::MarkovChain, s0, t0::DateTime, nstep::Int64; ny::Int64=1, ns::Int64=1)
    # Parameters
    n = size(mc_wk.states[1],2)
    n_state = size(mc_wk.transition_matrices[1],1)
    t = t0:Hour(1):t0+Hour(nstep-1)

    # Pre-allocate
    ω_generated = [(t = repeat(t, 1, ny, ns), power = zeros(nstep, ny, ns)) for _ in 1:n]
    probabilities = zeros(ny, ns)

    # Initialization
    for j in 1:n
        ω_generated[j].power[1,:,:] .= s0[j]
    end

    # Initialize state with the closest value
    isweekend(t0) ? mc = mc_wkd : mc = mc_wk
    Δ_s0 = [s0 .- mc.states[Dates.month(t0)][Dates.hour(t0)+1, :, :][:,state] for state in 1:n_state]
    idx_0 = findmin(norm.(Δ_s0))[2]

    # Simulate
    for s in 1:ns, y in 1:ny, k in 2:nstep
        # Retrieve the current hour and month - hour(t[1]) = 0...
        h, m = Dates.hour(t[k]) + 1, Dates.month(t[k])
        # Test to chose the appropriate markov chain
        isweekend(t[k]) ? mc = mc_wkd : mc = mc_wk
        # Build categorical distribution with transition matrices
        distributions = [Categorical(mc.transition_matrices[m][state, :, max(1, h - 1)]) for state in 1:n_state]
        # Compute state index from the probability matrix
        idx_1 = rand(distributions[idx_0])
        # Retrieve the associated values
        for j in 1:n
            ω_generated[j].power[k,y,s] = mc.states[m][h, j, idx_1]
        end
        # Update idx
        idx_0 = idx_1
    end

    return ω_generated
end
# Generate scenario from a single markov chain
function generate(generator::MarkovGeneration, s0, t0::DateTime, nstep::Int64; ny::Int64=1, ns::Int64=1)
    # Parameters
    n = size(generator.mc.states[1],2)
    n_state = size(generator.mc.transition_matrices[1],1)
    t = t0:Hour(1):t0+Hour(nstep-1)

    # Pre-allocate
    ω_generated = [(t = repeat(t, 1, ny, ns), power = zeros(nstep, ny, ns)) for _ in 1:n]
    probabilities = zeros(ny, ns)

    # Initialization
    for j in 1:n
        ω_generated[j].power[1,:,:] .= s0[j]
    end

    # Initialize state with the closest value
    Δ_s0 = [s0 .- generator.mc.states[Dates.month(t0)][Dates.hour(t0)+1, :, :][:,state] for state in 1:n_state]
    idx_0 = findmin(norm.(Δ_s0))[2]

    # Simulate
    for s in 1:ns, y in 1:ny, k in 2:nstep
        # Retrieve the current hour and month - hour(t[1]) = 0...
        h, m = Dates.hour(t[k]) + 1, Dates.month(t[k])
        # Build categorical distribution with transition matrices
        distributions = [Categorical(generator.mc.transition_matrices[m][state, :, max(1, h - 1)]) for state in 1:n_state]
        # Compute state index from the probability matrix
        idx_1 = rand(distributions[idx_0])
        # Retrieve the associated values
        for j in 1:n
            ω_generated[j].power[k,y,s] = generator.mc.states[m][h, j, idx_1]
        end
        # Update idx
        idx_0 = idx_1
    end

    return ω_generated, probabilities
end

# Generate a yearly profile from typical days clustered data
generate(data_td::Array{Float64,2}, assignments::Array{Int64,1}) = reshape(data_td[:, assignments, :], size(data_td, 1) * size(assignments, 1))
