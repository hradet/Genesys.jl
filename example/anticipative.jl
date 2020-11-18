# Load packages
using Genesys, CSV, DataFrames, JLD, Dates, Seaborn
pygui(true)

# Parameters
const nh, ny, ns = 8760, 5, 1

# Load data
data = load(joinpath("data","scenarios_ausgrid.jld"))

# Initialize scenarios
ω_simu = Scenarios(data["ω_simu"], 1:nh, 1:ny, 1:ns)

# Initialize DES
des = DistributedEnergySystem(ld_E = Load(),
                              pv = Source(),
                              liion = Liion(),
                              grid = Grid(),
                              parameters = Genesys.GlobalParameters(nh, ny, ns, τ_share = 0.8))

# Offline computations
controller, designer = offline_optimization!(des,
                                            AnticipativeTwoStage(),
                                            ω_simu)

# Simulate
@elapsed simulate!(des, controller, designer, ω_simu)
