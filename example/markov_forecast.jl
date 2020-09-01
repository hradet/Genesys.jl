using Genesys, Statistics, CSV, Seaborn
pygui(true)

# Parameters
n_markov = 20 # number of markov state
nω = 100 # number of scenario to simulate
nstep = 24 # length of the scenarios
h, y, s = 1, 1, 1 # year studied

# GUI function
include("load_GUI.jl")

# GUI loading
outputGUI = loadGUI([], [], 1, 1, 1)

# Initialization
ld, pv, liion, h2tank, elyz, fc, tes, heater, controller, designer, grid, ω_optim, ω_simu = Genesys.initialization(outputGUI)

# Clustering data
pv_month = Genesys.clustering_month(ω_optim.pv_E, ω_optim.timestamp)
ld_wk = Genesys.clustering_month_week(ω_optim.ld_E, ω_optim.timestamp)
ld_wkd = Genesys.clustering_month_weekend(ω_optim.ld_E, ω_optim.timestamp)

# Compute markov chain
mc_pv = Genesys.compute_markovchain(pv_month, n_markov)
mc_ld_wk = Genesys.compute_markovchain(ld_wk, n_markov)
mc_ld_wkd = Genesys.compute_markovchain(ld_wkd, n_markov)

# Compute scenario forecast from MC
ω_pv = zeros(nstep, nω)
ω_ld = zeros(nstep, nω)
for ω in 1:nω
    ω_pv[:,ω] = Genesys.compute_scenario(mc_pv, ω_optim.pv_E[h,y,s], ω_optim.timestamp[h], y, nstep-1)
    ω_ld[:,ω] = Genesys.compute_scenario(mc_ld_wk, mc_ld_wkd, ω_optim.ld_E[h,y,s], ω_optim.timestamp[h], y, nstep-1)
end

# Plots
plot(mean(ω_pv,dims=2), color="blue")
plot(mean(reshape(ω_optim.pv_E[1:24*31,1,1],24,:),dims=2), color="darkred")

# plot(mean(ω_ld,dims=2), color="blue")
# plot(ω_optim.ld_E[1:24,1,1], color="darkred")
