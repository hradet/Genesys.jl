module Genesys

# Optimisation
using JuMP, CPLEX, Metaheuristics, SDDP
# Math
using Statistics, StatsBase, MultivariateStats, Clustering, Distributions, Distances, LinearAlgebra
# using UMAP
# Others
using Seaborn, ProgressMeter, Dates, Distributed, SharedArrays, CSV, DataFrames, JLD
# Assets
include(joinpath("assets","microgrid.jl"))
include(joinpath("assets","carriers.jl"))
include(joinpath("assets","liion.jl"))
include(joinpath("assets","tes.jl"))
include(joinpath("assets","h2tank.jl"))
include(joinpath("assets","electrolyzer.jl"))
include(joinpath("assets","fuelcell.jl"))
include(joinpath("assets","heater.jl"))
include(joinpath("assets","grid.jl"))
include(joinpath("assets","solar.jl"))
include(joinpath("assets","demand.jl"))
export Microgrid, Demand, Solar, Liion, ThermalStorage, H2Tank, FuelCell, Electrolyzer, Heater, Grid, GlobalParameters
export Electricity, Heat, Hydrogen
export add!
# # Scenarios
include(joinpath("scenarios","scenarios.jl"))
include(joinpath("scenarios","reduction.jl"))
include(joinpath("scenarios","generation.jl"))
include(joinpath("scenarios","utils.jl"))
export Scenarios
export ManualReducer, SAAReducer, MeanValueReducer, FeatureBasedReducer
export UnitRangeTransform, ZScoreTransform
export PCAReduction, StatsReduction
export KmedoidsClustering
export MarkovGenerator, AnticipativeGenerator
export reduce, generate
# # Optimization utils
include(joinpath("optimization","utils.jl"))
export Expectation, CVaR, WorstCase
# # Operation optimization
include(joinpath("optimization","controller","dummy.jl"))
include(joinpath("optimization","controller","anticipative.jl"))
include(joinpath("optimization","controller","rb.jl"))
# include(joinpath("optimization","controller","olfc.jl"))
# include(joinpath("optimization","controller","sddp.jl"))
export Dummy, RBC, Anticipative#, OLFC, SDDPC
export RBCOptions, AnticipativeOptions#, OLFCOptions, SDDPCOptions
export initialize_controller!
# # Investment optimization
include(joinpath("optimization","designer","manual.jl"))
include(joinpath("optimization","designer","milp.jl"))
include(joinpath("optimization","designer","metaheuristic.jl"))
export Manual, Metaheuristic, MILP
export MetaheuristicOptions, MILPOptions
export initialize_designer!
# # Anticipative optimization
# include(joinpath("optimization","anticipative","multistages.jl"))
# include(joinpath("optimization","anticipative","twostage.jl"))
# export AnticipativeMultiStages, AnticipativeTwoStage
# export offline_optimization!
# # Simulation
include(joinpath("simulation","informations.jl"))
include(joinpath("simulation","dynamics.jl"))
include(joinpath("simulation","power_balances.jl"))
include(joinpath("simulation","simulations.jl"))
export simulate!
# # Utils
include(joinpath("utils","metrics.jl"))
include(joinpath("utils","plots.jl"))
# include(joinpath("utils","saves.jl"))
export Metrics
export plot_operation, plot_statistics

end
