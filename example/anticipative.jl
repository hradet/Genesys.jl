# Load packages
using Genesys, CSV, DataFrames, JLD, Dates, Seaborn
pygui(true)

# Parameters
const nh, ny, ns = 8760, 5, 1

# Load data
data = load(joinpath("data","input_data_stochastic.jld"))

# Initialize scenarios
ω_simu = Scenarios(data["ω_simu"], nh, ny, ns)

# Initialize DES
DES = DistributedEnergySystem(ld_E = Load(),
                              pv = Source(),
                              liion = Liion(),
                              grid = Grid(),
                              parameters = Genesys.GlobalParameters(nh, ny, ns, τ_share = 0.8))

# Offline computations
controller, designer = offline_optimization!(DES,
                                            AnticipativeEAC(options = Genesys.EACStochOptions(range_y = 1:ny)),
                                            ω_simu)

# Simulate
@elapsed simulate!(DES, controller, designer, ω_simu)
