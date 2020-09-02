module Genesys

# Optimisation
using JuMP, CPLEX
# Math
using Statistics, Clustering, Distributions
# Others
using Seaborn, ProgressMeter, Dates, Distributed, SharedArrays
# Components
include(joinpath("simulation","components","liion.jl"))
include(joinpath("simulation","components","tes.jl"))
include(joinpath("simulation","components","h2tank.jl"))
include(joinpath("simulation","components","electrolyzer.jl"))
include(joinpath("simulation","components","fuel_cell.jl"))
include(joinpath("simulation","components","heater.jl"))
include(joinpath("simulation","components","grid.jl"))
include(joinpath("simulation","components","sources.jl"))
include(joinpath("simulation","components","loads.jl"))
include(joinpath("simulation","components","abstract.jl"))
# Simulation
include(joinpath("simulation","informations.jl"))
include(joinpath("simulation","dynamics.jl"))
include(joinpath("simulation","power_balances.jl"))
include(joinpath("simulation","simulations.jl"))
# Scenarios
include(joinpath("scenarios","struct.jl"))
include(joinpath("scenarios","utils.jl"))
include(joinpath("scenarios","markov_functions.jl"))
include(joinpath("scenarios","time_blocks_functions.jl"))
include(joinpath("scenarios","typical_days_functions.jl"))
include(joinpath("scenarios","initialize_scenarios.jl"))
# Anticipative optimization
include(joinpath("optimization","anticipative","anticipative_controller.jl"))
include(joinpath("optimization","anticipative","anticipative_multistage.jl"))
include(joinpath("optimization","anticipative","anticipative_onestage.jl"))
include(joinpath("optimization","anticipative","anticipative_onestage_online_update.jl"))
include(joinpath("optimization","anticipative","utils.jl"))
# Investment optimization
include(joinpath("optimization","investment","rule_based.jl"))
include(joinpath("optimization","investment","eac.jl"))
include(joinpath("optimization","investment","eac_stoch.jl"))
# Operation optimization
include(joinpath("optimization","operation","rule_based.jl"))
include(joinpath("optimization","operation","mpc.jl"))
# Post-processing
include(joinpath("postprocessing","economics.jl"))
include(joinpath("postprocessing","plots.jl"))
include(joinpath("postprocessing","save_functions.jl"))
# Utils
include(joinpath("utils","initialization.jl"))

end
