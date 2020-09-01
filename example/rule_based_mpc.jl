# Genesys
using Genesys
# Lecture CSV
using JLD, Dates
# Plot
using Seaborn
pygui(true)

# GUI function
include("load_GUI.jl")

# Global constant
controller_flag, designer_flag = "MPCController" , "RuleBasedDesigner" # design and control methods
τ_energy_grid, τ_power_grid, τ_cost_grid= 1., 0., 1. # grid constraints

# GUI loading
outputGUI = loadGUI(controller_flag, designer_flag, τ_cost_grid, τ_power_grid, τ_energy_grid)

# Initialization
ld, pv, liion, h2tank, elyz, fc, tes, heater, controller, designer, grid, ω_optim, ω_simu = Genesys.initialization(outputGUI)

# Initialize designer
Genesys.initialize_designer(ld, pv, liion, designer, grid, ω_optim, outputGUI["parameters"])

# Initialize controller
Genesys.initialize_controller(ld, pv, liion, controller, grid, ω_optim, outputGUI["parameters"])

# Simulate
timer_online = @elapsed costs = Genesys.simulate(ld, pv, liion, controller, designer, grid, ω_optim, ω_simu, outputGUI["parameters"])

# Plot
Genesys.plot_operation(ld, pv, liion, grid, outputGUI["parameters"])
Genesys.plot_investment(designer, outputGUI["parameters"])
Genesys.plot_soh(liion)
