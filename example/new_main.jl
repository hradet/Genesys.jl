# Load packages
using Genesys, JLD, Dates, Seaborn
pygui(true)

# Parameters
const nh, ny, ns = 8760, 2, 1

# Load data
data = load(joinpath("data","ausgrid_5_twostage.jld"))

# Initialize scenarios
ω_optim, ω_simu = Scenarios(data["ω_optim"], 1:nh, 1:ny, 1:ns), Scenarios(data["ω_simu"],  1:nh, 1:ny, 1:ns)

# Create microgrid
microgrid = Microgrid(parameters = GlobalParameters(nh, ny, ns, renewable_share = 0.9999))

# Build the microgrid
add!(microgrid, Liion(), Solar(), Demand(carrier = Electricity()), Grid(carrier = Electricity()))

# Initialize controller

# Initialize designer

# Simulate
@elapsed simulate!(microgrid, controller, designer, ω_simu, options = Genesys.Options(mode="serial"))

# Compute the metrics
metrics = Metrics(microgrid, designer)
