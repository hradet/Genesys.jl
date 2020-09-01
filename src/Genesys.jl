module Genesys

# Optimisation
using JuMP, CPLEX
# Math
using Statistics, Clustering, Distributions
# Others
using Seaborn, ProgressMeter, Dates, Distributed, SharedArrays
#                                Simulation
#_______________________________________________________________________________
# Models
# Storages
include(joinpath("simulation","models","liion.jl"))
include(joinpath("simulation","models","tes.jl"))
include(joinpath("simulation","models","h2tank.jl"))
# Energy converters
include(joinpath("simulation","models","electrolyzer.jl"))
include(joinpath("simulation","models","fuel_cell.jl"))
include(joinpath("simulation","models","heater.jl"))
# Grid
include(joinpath("simulation","models","grid.jl"))
# Sources
include(joinpath("simulation","models","sources.jl"))
# Loads
include(joinpath("simulation","models","loads.jl"))
# Others
include(joinpath("simulation","struct.jl"))

# Multi-energy DES
include(joinpath("simulation","multi_energy","compute_economics.jl"))
include(joinpath("simulation","multi_energy","information_functions.jl"))
include(joinpath("simulation","multi_energy","dynamic_functions.jl"))
include(joinpath("simulation","multi_energy","checking_functions.jl"))
include(joinpath("simulation","multi_energy","simulate_functions.jl"))

# Simple DES
include(joinpath("simulation","simple","compute_economics.jl"))
include(joinpath("simulation","simple","information_functions.jl"))
include(joinpath("simulation","simple","dynamic_functions.jl"))
include(joinpath("simulation","simple","checking_functions.jl"))
include(joinpath("simulation","simple","simulate_functions.jl"))

#                                Scenarios
#_______________________________________________________________________________
include(joinpath("scenarios","struct.jl"))
include(joinpath("scenarios","utils.jl"))
include(joinpath("scenarios","markov_functions.jl"))
include(joinpath("scenarios","time_blocks_functions.jl"))
include(joinpath("scenarios","typical_days_functions.jl"))
include(joinpath("scenarios","initialize_scenarios.jl"))

#                                Optimization
#_______________________________________________________________________________
include(joinpath("optimization","utils.jl"))

# Operation
include(joinpath("optimization","operation","struct.jl"))
# Multi-energy DES
include(joinpath("optimization","operation","multi_energy","rule_based.jl"))
include(joinpath("optimization","operation","multi_energy","anticipative.jl"))
# Simple DES
include(joinpath("optimization","operation","simple","rule_based.jl"))
include(joinpath("optimization","operation","simple","anticipative.jl"))
include(joinpath("optimization","operation","simple","mpc.jl"))

# Investment
include(joinpath("optimization","investment","struct.jl"))
# Multi-energy DES
include(joinpath("optimization","investment","multi_energy","rule_based.jl"))
include(joinpath("optimization","investment","multi_energy","anticipative_multistage.jl"))
include(joinpath("optimization","investment","multi_energy","anticipative_onestage.jl"))
# Simple DES
include(joinpath("optimization","investment","simple","rule_based.jl"))
include(joinpath("optimization","investment","simple","anticipative_multistage.jl"))
include(joinpath("optimization","investment","simple","anticipative_onestage.jl"))
include(joinpath("optimization","investment","simple","anticipative_onestage_online_update.jl"))

#                                Utils
#_______________________________________________________________________________
include(joinpath("utils","initialization.jl"))
include(joinpath("utils","plot_functions.jl"))
include(joinpath("utils","save_functions.jl"))

end
