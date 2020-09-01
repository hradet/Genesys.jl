using Genesys, Statistics, Seaborn
pygui(true)

# Parameters
horizon = 100 # scenario length
states = 1:10 # markov states
nstates = length(s) # number of markov-states

# Create scenario randomly from the states
scenario = rand(states, horizon)

# Compute matrix transition
transition_matrix = Genesys.compute_transition_from_sequence(scenario, nstates)

# Create markov chain to compute a scenario
mc = Genesys.MarkovChain(states, transition_matrix)

# Compute scenario
scenario_markov = Genesys.compute_scenario(mc, states[1], horizon)
