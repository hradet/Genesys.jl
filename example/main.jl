# Genesys
using Genesys
# Lecture CSV
using CSV, DataFrames, JLD, Dates
# Plot
using Seaborn
pygui(true)

# Parameters
const nh, ny, ns = 8760, 20, 1

# Load data
data = load(joinpath("data","input_data_stochastic.jld"))

# Initialize scenarios
ω_optim, ω_simu = Scenarios(data["ω_optim"], nh, ny, ns), Scenarios(data["ω_simu"],  nh, ny, ns)

# Initialize DES
DES = DistributedEnergySystem(ld_E = Load(),
                              ld_H = Load(),
                              pv = Source(),
                              liion = Liion(),
                              tes = ThermalSto(),
                              h2tank = H2Tank(),
                              elyz = Electrolyzer(),
                              fc = FuelCell(),
                              heater = Heater(),
                              grid = Grid(),
                              parameters = Genesys.GlobalParameters(nh = nh, ny = ny, ns = ns, τ_share = 0.9))

# Initialize controller
controller = initialize_controller!(DES,
                                    RuleBasedController(),
                                    ω_optim)

# Initialize designer
designer = initialize_designer!(DES,
                                EACDesigner(),
                                ω_optim)

# Simulate
@elapsed simulate!(DES, controller, designer, ω_simu,
                          options = Genesys.Options(mode="serial"))

# Postprocessing
costs = Genesys.compute_economics(DES, designer)
tech = Genesys.compute_tech_indicators(DES)
Genesys.plot_operation(DES)
Genesys.plot_investment(DES, designer)
Genesys.plot_soh(DES)
Genesys.plot_costs(costs)
