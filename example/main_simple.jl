# Genesys
using Genesys
# Lecture CSV
using CSV, DataFrames, JLD, Dates
# Plot
using Seaborn
pygui(true)

# GUI loading
outputGUI = Genesys.loadGUI("RuleBasedController", "MetaHeuristicDesigner", τ_energy = 0.8)

# Initialization
ld, pv, liion, _, _, _, _, _, controller, designer, grid, ω_optim, ω_simu = Genesys.initialization(outputGUI)

# Initialize controller
Genesys.initialize_controller(ld, pv, liion, controller, grid, ω_optim, outputGUI["parameters"])

# Initialize designer
timer = @elapsed Genesys.initialize_designer(ld, pv, liion, controller, designer, grid, ω_optim, outputGUI["parameters"])

# Simulate
timer_online = @elapsed Genesys.simulate(ld, pv, liion, controller, designer, grid, ω_optim, ω_simu, outputGUI["parameters"], mode = "serial")

# Postprocessing
costs = Genesys.compute_economics(ld, pv, liion, designer, grid, outputGUI["parameters"])
tech = Genesys.compute_tech_indicators(ld, grid)
Genesys.plot_operation(ld, pv, liion, grid, outputGUI["parameters"])
Genesys.plot_investment(designer, outputGUI["parameters"])
Genesys.plot_soh(liion)
Genesys.plot_npv(costs)
