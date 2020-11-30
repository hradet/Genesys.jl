module Genesys

# Optimisation
using JuMP, CPLEX, Metaheuristics
# Math
using Statistics, Clustering, Distributions, Distances, LinearAlgebra
# Others
using Seaborn, ProgressMeter, Dates, Distributed, SharedArrays, CSV, DataFrames, JLD
# Components
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
export ManualReducer, SAAReducer, ExpectedValueReducer, KmeansReducer, KmedoidsReducer
export MarkovGenerator
export reduce, generate
# Operation optimization
include(joinpath("optimization","controller","abstract.jl"))
include(joinpath("optimization","controller","dummy.jl"))
include(joinpath("optimization","controller","anticipative.jl"))
include(joinpath("optimization","controller","rb.jl"))
include(joinpath("optimization","controller","olfc.jl"))
export DummyCcontroller, RBC, OLFC, Anticipative
export initialize_controller!
# Investment optimization
include(joinpath("optimization","designer","abstract.jl"))
include(joinpath("optimization","designer","dummy.jl"))
include(joinpath("optimization","designer","milp.jl"))
include(joinpath("optimization","designer","metaheuristic.jl"))
export DummyDesigner, MILP, Metaheuristic
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
include(joinpath("utils","GUI.jl"))
export Metrics
export plot_operation, plot_investment, plot_soh, plot_costs, plot_statistics

end
