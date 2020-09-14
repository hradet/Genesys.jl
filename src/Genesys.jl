module Genesys

# Optimisation
using JuMP, CPLEX, Metaheuristics
# Math
using Statistics, Clustering, Distributions
# Others
using Seaborn, ProgressMeter, Dates, Distributed, SharedArrays, CSV, DataFrames, JLD
# Components
include(joinpath("simulation","items","liion.jl"))
include(joinpath("simulation","items","tes.jl"))
include(joinpath("simulation","items","h2tank.jl"))
include(joinpath("simulation","items","electrolyzer.jl"))
include(joinpath("simulation","items","fuel_cell.jl"))
include(joinpath("simulation","items","heater.jl"))
include(joinpath("simulation","items","grid.jl"))
include(joinpath("simulation","items","sources.jl"))
include(joinpath("simulation","items","loads.jl"))
include(joinpath("simulation","items","abstract.jl"))
include(joinpath("simulation","items","DES.jl"))
# Simulation
include(joinpath("simulation","informations.jl"))
include(joinpath("simulation","dynamics.jl"))
include(joinpath("simulation","power_balances.jl"))
include(joinpath("simulation","simulations.jl"))
# Scenarios
include(joinpath("scenarios","scenario_reductions.jl"))
include(joinpath("scenarios","clustering.jl"))
include(joinpath("scenarios","markov.jl"))
include(joinpath("scenarios","utils.jl"))
# Anticipative optimization
include(joinpath("optimization","anticipative","anticipative_controller.jl"))
include(joinpath("optimization","anticipative","anticipative_multistage.jl"))
include(joinpath("optimization","anticipative","anticipative_onestage.jl"))
include(joinpath("optimization","anticipative","anticipative_onestage_online_update.jl"))
# Investment optimization
include(joinpath("optimization","designer","dummy.jl"))
include(joinpath("optimization","designer","rule_based.jl"))
include(joinpath("optimization","designer","eac.jl"))
include(joinpath("optimization","designer","eac_td.jl"))
include(joinpath("optimization","designer","eac_stoch.jl"))
include(joinpath("optimization","designer","eac_stoch_td.jl"))
include(joinpath("optimization","designer","metaheuristic.jl"))
# Operation optimization
include(joinpath("optimization","controller","dummy.jl"))
include(joinpath("optimization","controller","rule_based.jl"))
include(joinpath("optimization","controller","mpc.jl"))
# Post-processing
include(joinpath("postprocessing","indicators.jl"))
include(joinpath("postprocessing","plots.jl"))
include(joinpath("postprocessing","saves.jl"))
# Utils
include(joinpath("utils","GUI.jl"))

end
