# Load packages
using Genesys, JLD, Dates, Seaborn
pygui(true)

# Parameters
const nh, ny, ns = 8760, 2, 10

# Load input data
data = load(joinpath("data","ausgrid_5_twostage.jld"))

# Create microgrid
microgrid = Microgrid(parameters = GlobalParameters(nh, ny, ns, renewable_share = 0.9999))

# Build the microgrid
add!(microgrid, Demand(carrier = Electricity()), Demand(carrier = Heat()),
                Solar(),
                Liion(), ThermalStorage(), H2Tank(),
                Heater(), Electrolyzer(), FuelCell(),
                Grid(carrier = Electricity()))

# Initialize scenarios
ω_optim, ω_simu = Scenarios(data["ω_optim"], microgrid), Scenarios(data["ω_simu"],  microgrid)

# Initialize controller
controller = initialize_controller!(microgrid, RBC(options = RBCOptions(policy_selection = 1)), ω_optim)

# Initialize designer
designer = initialize_designer!(microgrid, Manual(generations = [112.], storages = [149., 585., 1597.], converters = [30., 1.5, 3.3]), ω_optim)

designer = initialize_designer!(microgrid, Metaheuristic(options = MetaheuristicOptions(reducer = ManualReducer(), iterations = 1)), ω_optim)

# Simulate
@elapsed simulate!(microgrid, controller, designer, ω_simu, options = Genesys.Options(mode="serial"))

# Compute the metrics
metrics = Metrics(microgrid, designer)

# Plots
plot_operation(microgrid)
plot_statistics(metrics)

heat = sum(Genesys.power_balance(1:nh, 1:ny, 1:ns, microgrid, typeof(Heat())), dims=1)
