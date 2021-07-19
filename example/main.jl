# Load packages
using Genesys, JLD, Dates, Seaborn
pygui(true)

# Parameters
const nh, ny, ns = 8760, 2, 1

# Load data
data = load(joinpath("data","ausgrid_5_twostage.jld"))

# Initialize scenarios
ω_optim, ω_simu = Scenarios(data["ω_optim"], 1:nh, 1:ny, 1:ns), Scenarios(data["ω_simu"],  1:nh, 1:ny, 1:ns)

# Initialize DES
des = DistributedEnergySystem(ld_E = Load(),
                              pv = Source(),
                              liion = Liion(),
                              grid = Grid(),
                              parameters = Genesys.GlobalParameters(nh, ny, ns, renewable_share = 0.8))

# Initialize controller
controller = initialize_controller!(des,
                                    RBC(),
                                    ω_optim)

# Initialize designer
designer = initialize_designer!(des,
                                MILP(options = MILPOptions(reducer = ManualReducer())),
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
