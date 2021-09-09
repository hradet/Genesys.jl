# Load packages
using Genesys, JLD, Dates, Seaborn, BenchmarkTools
pygui(true)

# Parameters
const nh, ny, ns = 8760, 2, 1

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

designer = initialize_designer!(microgrid, Metaheuristic(options = MetaheuristicOptions(controller = RBC(options = RBCOptions(policy_selection = 3)), reducer = ManualReducer(y = 2:2, s = 179:179), iterations = 20, multithreads = false)), ω_optim)

designer = initialize_designer!(microgrid, MILP(options = MILPOptions(reducer = ManualReducer(y=1:1, s=1:1))), ω_optim)

controller = initialize_controller!(microgrid, Anticipative(generations = [a[1] for a in designer.decisions.generations], storages = [a[1] for a in designer.decisions.storages], converters = [a[1] for a in designer.decisions.converters]), ω_optim)

options = OLFCOptions(nscenarios = 5, reducer = ManualReducer(y=1:1, s=1:1))
controller = initialize_controller!(microgrid, OLFC(options = options, generations = [a[1] for a in designer.decisions.generations], storages = [a[1] for a in designer.decisions.storages], converters = [a[1] for a in designer.decisions.converters]), ω_optim)

# Simulate
@time simulate!(microgrid, controller, designer, ω_optim, options = Genesys.Options(mode="serial"))

# Compute the metrics
metrics = Metrics(microgrid, designer)

# Plots
plot_operation(microgrid)
plot_metrics(metrics)

a,b,c,d = Genesys.compute_forecast(1, 2, 1, microgrid, controller)
