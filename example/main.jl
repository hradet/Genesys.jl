# Load packages
using Genesys, JLD, Dates, Seaborn
pygui(true)

# Parameters
const nh, ny, ns = 8760, 2, 100

# Load input data
data = load(joinpath("data","ausgrid_5_twostage.jld"))

# Create microgrid
microgrid = Microgrid(parameters = GlobalParameters(nh, ny, ns, renewable_share = 1.))

# Build the microgrid
add!(microgrid, Demand(carrier = Electricity()), Demand(carrier = Heat()),
                Solar(),
                Liion(), ThermalStorage(), H2Tank(),
                Heater(), Electrolyzer(), FuelCell(),
                Grid(carrier = Electricity()))

# Initialize scenarios
ω_d, ω_a = Scenarios(data["ω_optim"], microgrid), Scenarios(data["ω_simu"], microgrid)

# Initialize controller
controller = initialize_controller!(microgrid, RBC(), ω_d)

# Initialize designer
designer = initialize_designer!(microgrid, MILP(options = MILPOptions(reducer = ManualReducer(y=1:1, s=1:1))), ω_d)

# Simulate
simulate!(microgrid, controller, designer, ω_a, options = Genesys.Options(mode = "multithreads"))

# Compute the metrics
metrics = Metrics(microgrid, designer)

# Plots
plot_operation(microgrid)
plot_metrics(metrics)
