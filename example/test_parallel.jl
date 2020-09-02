# Lecture CSV
using CSV, DataFrames, JLD, Dates
# Plot
using Seaborn
pygui(true)

using Distributed
addprocs(2)
@everywhere using Genesys

# GUI function
include("load_GUI.jl")

# GUI loading
outputGUI = loadGUI("RuleBasedController", "RuleBasedDesigner", 1., 0., 1.)

# Initialization
ld, pv, liion, h2tank, elyz, fc, tes, heater, controller, designer, grid, ω_optim, ω_simu = Genesys.initialization(outputGUI)

# Initialize designer
Genesys.initialize_designer(ld, pv, liion, designer, grid, ω_optim, outputGUI["parameters"])

# Initialize controller
Genesys.initialize_controller(ld, pv, liion, controller, grid, ω_optim, outputGUI["parameters"])

# Simulate serial
timer_online_s = @elapsed costs_s = Genesys.simulate(ld, pv, liion, controller, designer, grid, ω_optim, ω_simu, outputGUI["parameters"], mode="serial")

# Simulate parallel
timer_online_p1 = @elapsed costs_p1 = Genesys.simulate(ld, pv, liion, controller, designer, grid, ω_optim, ω_simu, outputGUI["parameters"], mode="multicores")

# Simulate parallel distributed
timer_online_p2 = @elapsed costs_p2 = Genesys.simulate(ld, pv, liion, controller, designer, grid, ω_optim, ω_simu, outputGUI["parameters"], mode="distributed")

# Simulate multi-threads
timer_online_p3 = @elapsed costs_p3 = Genesys.simulate(ld, pv, liion, controller, designer, grid, ω_optim, ω_simu, outputGUI["parameters"], mode="multithreads")

# Plot
Genesys.plot_operation(ld, pv, liion, grid, outputGUI["parameters"])
Genesys.plot_investment(designer, outputGUI["parameters"])
Genesys.plot_soh(liion)
Genesys.plot_npv(costs)
