# Load packages
using Genesys, CSV, DataFrames, JLD, Dates, Seaborn
pygui(true)

# Parameters
const nh, ny, ns = 8760, 20, 10

# Load data
data = load(joinpath("data","scenarios_ausgrid.jld"))

# Initialize scenarios
ω_optim, ω_simu = Scenarios(data["ω_optim"], 1:nh, 1:ny, 1:ns), Scenarios(data["ω_simu"],  1:nh, 1:ny, 1:ns)

# Initialize DES
des = DistributedEnergySystem(ld_E = Load(),
                              pv = Source(),
                              liion = Liion(),
                              grid = Grid(),
                              parameters = Genesys.GlobalParameters(nh, ny, ns, τ_share = 0.5))

# Initialize controller
controller = initialize_controller!(des,
                                    RBC(),
                                    ω_optim)

# Initialize designer
designer = initialize_designer!(des,
                                MILP(),
                                ω_optim)

# Simulate
@elapsed simulate!(des, controller, designer, ω_simu,
                          options = Genesys.Options(mode="serial"))

# Compute the metrics
metrics = Metrics(des, designer)

# Plot functions
plot_operation(des)
plot_investment(des, designer)
plot_soh(des)
plot_costs(metrics.costs)
plot_statistics(metrics)
