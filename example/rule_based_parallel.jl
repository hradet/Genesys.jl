# Lecture CSV
using CSV, DataFrames, JLD, Dates
# Plot
using Seaborn
pygui(true)

#TODO: Add to Genesys.jl as a function
#TODO: Load Genesys as a local package
using Distributed
addprocs(2)
@everywhere begin
 push!(LOAD_PATH, "C:\\Users\\radet\\SynologyDrive\\01 THESE\\01 ESPACE DEVELOPPEMENT\\Julia\\Genesys.jl\\src")
 @eval using Genesys
end

# GUI function
include("load_GUI.jl")
# Global constant
controller_flag, designer_flag = "RuleBasedController" , "RuleBasedDesigner" # design and control methods
τ_energy_grid, τ_power_grid, τ_cost_grid= 1., 0., 1. # grid constraints


# GUI loading
outputGUI = loadGUI(controller_flag, designer_flag, τ_cost_grid, τ_power_grid, τ_energy_grid)

# Initialization
ld, pv, liion, h2tank, elyz, fc, tes, heater, controller, designer, grid, ω_optim, ω_simu = Genesys.initialization(outputGUI)

# Initialize designer
Genesys.initialize_designer(ld, pv, liion, designer, grid, ω_optim, outputGUI["parameters"])

# Initialize controller
Genesys.initialize_controller(ld, pv, liion, controller, grid, ω_optim, outputGUI["parameters"])

# Simulate serial
timer_online_s = @elapsed costs_s = Genesys.simulate(ld, pv, liion, controller, designer, grid, ω_optim, ω_simu, outputGUI["parameters"])

# Simulate parallel
timer_online_p1 = @elapsed costs_p1 = Genesys.simulate_parallel(ld, pv, liion, controller, designer, grid, ω_optim, ω_simu, outputGUI["parameters"])

# Simulate parallel distributed
timer_online_p2 = @elapsed costs_p2 = Genesys.simulate_parallel_distributed(ld, pv, liion, controller, designer, grid, ω_optim, ω_simu, outputGUI["parameters"])

# Simulate multi-threads
timer_online_p3 = @elapsed costs_p3 = Genesys.simulate_parallel_thread(ld, pv, liion, controller, designer, grid, ω_optim, ω_simu, outputGUI["parameters"])

# Plot
Genesys.plot_operation(ld, pv, liion, grid, outputGUI["parameters"])
Genesys.plot_investment(designer, outputGUI["parameters"])
Genesys.plot_soh(liion)
Genesys.plot_npv(costs)
