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
export DistributedEnergySystem, Source, Liion, ThermalSto, H2Tank, Electrolyzer, FuelCell, Heater, Grid, Load
# Simulation
include(joinpath("simulation","informations.jl"))
include(joinpath("simulation","dynamics.jl"))
include(joinpath("simulation","power_balances.jl"))
include(joinpath("simulation","simulations.jl"))
export simulate!
# Optimization
include(joinpath("optimization","utils.jl"))
# Investment optimization
include(joinpath("optimization","designer","dummy.jl"))
include(joinpath("optimization","designer","eac.jl"))
include(joinpath("optimization","designer","eac_stoch.jl"))
include(joinpath("optimization","designer","metaheuristic.jl"))
export DummyDesigner, EAC, EACStoch, Metaheuristic
export initialize_designer!
# Operation optimization
include(joinpath("optimization","controller","dummy.jl"))
include(joinpath("optimization","controller","rb.jl"))
include(joinpath("optimization","controller","mpc.jl"))
export DummyCcontroller, RBC, MPC
export initialize_controller!
# Anticipative optimization
include(joinpath("optimization","anticipative","controller.jl"))
include(joinpath("optimization","anticipative","multistages.jl"))
include(joinpath("optimization","anticipative","eac.jl"))
export Anticipative, AnticipativeMultiStages, AnticipativeEAC
export offline_optimization!
# Scenarios
include(joinpath("scenarios","scenario_reductions.jl"))
include(joinpath("scenarios","clustering.jl"))
include(joinpath("scenarios","markov.jl"))
include(joinpath("scenarios","utils.jl"))
export Scenarios
# Utils
include(joinpath("utils","metrics.jl"))
include(joinpath("utils","plots.jl"))
include(joinpath("utils","saves.jl"))
include(joinpath("utils","GUI.jl"))
export Metrics
export plot_operation, plot_investment, plot_soh, plot_costs, plot_statistics

end
