# Genesys
using Genesys
# Optimization
using JuMP
# Lecture CSV
using CSV, DataFrames, JLD, Dates
# Plot
using Seaborn
pygui(true)

#=
You can either use a controller and designer from the library,
or you can define your own controller and designer and modify the associated functions...
=#

# Define your own controller and designer
mutable struct foo <: Genesys.AbstractController
    u::NamedTuple
    foo() = new()
end
mutable struct bar <: Genesys.AbstractDesigner
    u::NamedTuple
    bar() = new()
end

# Define offline functions
function Genesys.initialize_controller!(des::DistributedEnergySystem, controller::foo, ω::Scenarios)
    # Preallocation
    Genesys.preallocate!(controller, des.parameters.nh, des.parameters.ny, des.parameters.ns)
    return controller
end

function Genesys.initialize_designer!(des::DistributedEnergySystem, designer::bar, ω::Scenarios)
    # Preallocation
    Genesys.preallocate!(designer, des.parameters.ny, des.parameters.ns)
    return designer
end

# Define online functions
function Genesys.compute_operation_decisions!(h::Int64, y::Int64, s::Int64, des::DistributedEnergySystem, controller::foo)
    return controller
end
function Genesys.compute_investment_decisions!(y::Int64, s::Int64, des::DistributedEnergySystem, designer::bar)
    return designer
end

#=
Let's simulate the microgrid with the dummies controller and designer...
=#

# Parameters
const nh, ny, ns = 8760, 20, 1

# Load data
data = load(joinpath("data","input_data_stochastic.jld"))

# Initialize scenarios
ω_optim, ω_simu = Scenarios(data["ω_optim"], nh, ny, ns), Scenarios(data["ω_simu"],  nh, ny, ns)

# Initialize DES
DES = DistributedEnergySystem(ld_E = Load(),
                              pv = Source(),
                              liion = Liion(),
                              grid = Grid(),
                              parameters = Genesys.GlobalParameters(nh = nh, ny = ny, ns = ns))

# Initialize controller
controller = initialize_controller!(DES,
                                    foo(),
                                    ω_optim)

# Initialize designer
designer = initialize_designer!(DES,
                                bar(),
                                ω_optim)

# Simulate
@elapsed simulate!(DES, controller, designer, ω_simu,
                          options = Genesys.Options(mode="serial"))
