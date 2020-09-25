# Load packages
using Genesys, CSV, DataFrames, JLD, Dates, Seaborn
pygui(true)

# Parameters
const nh, ny, ns = 8760, 20, 100

# Load data
data = load(joinpath("data","scenarios_ausgrid.jld"))

# Initialize scenarios
ω_optim, ω_simu = Scenarios(data["ω_optim"], nh, ny, ns), Scenarios(data["ω_simu"],  nh, ny, ns)

# Initialize DES
DES = DistributedEnergySystem(ld_E = Load(),
                              pv = Source(),
                              liion = Liion(),
                              grid = Grid(),
                              parameters = Genesys.GlobalParameters(nh = nh, ny = ny, ns = ns, τ_share = 0.8))

# Initialize controller
controller = initialize_controller!(DES,
                                    RBC(),
                                    ω_optim)

# Initialize designer
designer = initialize_designer!(DES,
                                EAC(),
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
