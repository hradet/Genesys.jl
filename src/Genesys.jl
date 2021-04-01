module Genesys

# Optimisation
using JuMP, CPLEX, Metaheuristics, SDDP
# Math
using Statistics, StatsBase, MultivariateStats, Clustering, Distributions, Distances, LinearAlgebra, HDBSCAN
# using UMAP
# Others
using Seaborn, ProgressMeter, Dates, Distributed, SharedArrays, CSV, DataFrames, JLD
# Assets
include(joinpath("assets","liion.jl"))
include(joinpath("assets","tes.jl"))
include(joinpath("assets","h2tank.jl"))
include(joinpath("assets","electrolyzer.jl"))
include(joinpath("assets","fuelcell.jl"))
include(joinpath("assets","heater.jl"))
include(joinpath("assets","grid.jl"))
include(joinpath("assets","sources.jl"))
include(joinpath("assets","loads.jl"))
include(joinpath("assets","des.jl"))
export DistributedEnergySystem, Source, Liion, ThermalSto, H2Tank, Electrolyzer, FuelCell, Heater, Grid, Load
# Scenarios
include(joinpath("scenarios","scenarios.jl"))
include(joinpath("scenarios","reduction.jl"))
include(joinpath("scenarios","generation.jl"))
include(joinpath("scenarios","utils.jl"))
export Scenarios
export ManualReducer, SAAReducer, MeanValueReducer, FeatureBasedReducer
export UnitRangeTransform, ZScoreTransform
export PCAReduction, StatsReduction
# export UMAPReduction
export KmedoidsClustering, HDBSCANClustering
export MarkovGenerator, AnticipativeGenerator
export reduce, generate
# Risk measures
include(joinpath("optimization","risk_measures.jl"))
export Expectation, CVaR, WorstCase
# Operation optimization
include(joinpath("optimization","controller","abstract.jl"))
include(joinpath("optimization","controller","dummy.jl"))
include(joinpath("optimization","controller","anticipative.jl"))
include(joinpath("optimization","controller","rb.jl"))
include(joinpath("optimization","controller","olfc.jl"))
include(joinpath("optimization","controller","sddp.jl"))
export Dummy, RBC, OLFC, Anticipative, SDDPC
export RBCOptions, OLFCOptions, AnticipativeOptions, SDDPCOptions
export initialize_controller!
# Investment optimization
include(joinpath("optimization","designer","abstract.jl"))
include(joinpath("optimization","designer","manual.jl"))
include(joinpath("optimization","designer","milp.jl"))
include(joinpath("optimization","designer","metaheuristic.jl"))
export Manual, MILP, Metaheuristic
export MILPOptions, MetaheuristicOptions
export initialize_designer!
# Anticipative optimization
include(joinpath("optimization","anticipative","multistages.jl"))
include(joinpath("optimization","anticipative","twostage.jl"))
export AnticipativeMultiStages, AnticipativeTwoStage
export offline_optimization!
# Simulation
include(joinpath("simulation","informations.jl"))
include(joinpath("simulation","dynamics.jl"))
include(joinpath("simulation","power_balances.jl"))
include(joinpath("simulation","simulations.jl"))
export simulate!
# Utils
include(joinpath("utils","metrics.jl"))
include(joinpath("utils","plots.jl"))
include(joinpath("utils","saves.jl"))
export Metrics
export plot_operation, plot_investment, plot_soh, plot_costs, plot_statistics

end
