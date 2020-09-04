# Genesys
using Genesys
# Lecture CSV
using CSV, DataFrames, JLD, Dates
# Plot
using Seaborn
pygui(true)

# GUI loading
outputGUI = Genesys.loadGUI("RuleBasedController", "RuleBasedDesigner")

# Initialization
ld, pv, liion, h2tank, elyz, fc, tes, heater, controller, designer, grid, ω_optim, ω_simu = Genesys.initialization(outputGUI)

# Initialize designer
Genesys.initialize_designer(ld, pv, liion, h2tank, elyz, fc, tes, heater, designer, grid, ω_optim, outputGUI["parameters"])

# Initialize controller
Genesys.initialize_controller(ld, pv, liion, h2tank, elyz, fc, tes, heater, controller, grid, ω_optim, outputGUI["parameters"])

# Simulate
timer_online = @elapsed Genesys.simulate(ld, pv, liion, h2tank, elyz, fc, tes, heater, controller, designer, grid, ω_optim, ω_simu, outputGUI["parameters"], mode = "serial")

# Postprocessing
costs = Genesys.compute_economics(ld, pv, liion, designer, grid, outputGUI["parameters"])
Genesys.plot_operation(ld, pv, liion, h2tank, elyz, fc, tes, heater, grid, outputGUI["parameters"])
Genesys.plot_investment(designer, outputGUI["parameters"])
Genesys.plot_soh(liion, elyz, fc)
Genesys.plot_npv(costs)
