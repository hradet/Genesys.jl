# Genesys
using Genesys
# Optimization
using JuMP
# Lecture CSV
using CSV, DataFrames, JLD2
# Plot
using Seaborn
pygui(true)

# GUI function
include("load_GUI.jl")

# Global constant
controller_flag, designer_flag = "AnticipativeController" , "AnticipativeOneStageDesigner" # design and control methods
τ_energy_grid, τ_power_grid, τ_cost_grid= 0.6, 0., 1. # grid constraints

# GUI loading
outputGUI = loadGUI(controller_flag, designer_flag, τ_cost_grid, τ_power_grid, τ_energy_grid)

# Initialization
ld, pv, liion, _, _, _, _, _, controller, designer, grid, ω_optim, ω_simu = Genesys.initialization(outputGUI)

# Offline optimization
timer_offline = @elapsed Genesys.offline_optimization(ld, pv, liion, controller, designer, grid, ω_optim, outputGUI["parameters"])

# Simulate
timer_online = @elapsed costs = Genesys.simulate(ld, pv, liion, controller, designer, grid, ω_optim, ω_simu, outputGUI["parameters"])

# Plot
Genesys.plot_operation(ld, pv, liion, grid, outputGUI["parameters"])
Genesys.plot_investment(designer, outputGUI["parameters"])
Genesys.plot_soh(liion)
