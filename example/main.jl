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
add!(microgrid, Demand(carrier = Electricity()), Demand(carrier = Heat()),
                Solar(),
                Liion(), ThermalStorage(), H2Tank(),
                Heater(), Electrolyzer(), FuelCell(),
                Grid(carrier = Electricity()))

# Initialize controller
controller = initialize_controller!(microgrid, RBC(options = RBCOptions(policy_selection = 1)), ω_optim)

# Initialize designer
designer = initialize_designer!(microgrid, Manual(generations = [112.], storages = [149., 585., 1597.], converters = [30., 1.5, 3.3]), ω_optim)

# Simulate
@elapsed simulate!(microgrid, controller, designer, ω_simu, options = Genesys.Options(mode="serial"))

# Compute the metrics
metrics = Metrics(microgrid, designer)

a = Genesys.baseline_cost(1:2, 1:1, microgrid)
