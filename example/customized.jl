# Packages
using Genesys, JLD, Dates

#=
You can either use a controller and designer from the library,
or you can define your own controller and designer and modify the associated functions...
=#

# Define your own controller and designer
mutable struct foo <: Genesys.AbstractController
    decisions::NamedTuple
    foo() = new()
end
mutable struct bar <: Genesys.AbstractDesigner
    decisions::NamedTuple
    bar() = new()
end

# Define offline functions
function Genesys.initialize_controller!(mg::Microgrid, controller::foo, ω::Scenarios)
    # Preallocate
    Genesys.preallocate!(mg, controller)
    return controller
end

function Genesys.initialize_designer!(mg::Microgrid, designer::bar, ω::Scenarios)
    # Preallocate
    Genesys.preallocate!(mg, designer)
    return designer
end

# Define online functions
function Genesys.compute_operation_decisions!(h::Int64, y::Int64, s::Int64, mg::Microgrid, controller::foo)
    return controller
end
function Genesys.compute_investment_decisions!(y::Int64, s::Int64, mg::Microgrid, designer::bar)
    return designer
end

#=
Let's simulate the microgrid with the dummies controller and designer...
=#

# Parameters
const nh, ny, ns = 8760, 2, 1

# Load input data
data = load(joinpath("data","ausgrid_5_twostage.jld"))

# Create microgrid
microgrid = Microgrid(parameters = GlobalParameters(nh, ny, ns))

# Build the microgrid
add!(microgrid, Demand(carrier = Electricity()), Solar(), Liion(), Grid(carrier = Electricity()))

# Initialize scenarios
ω_d, ω_a = Scenarios(data["ω_optim"], microgrid), Scenarios(data["ω_simu"],  microgrid)

# Initialize controller
controller = initialize_controller!(microgrid, foo(), ω_d)

# Initialize designer
designer = initialize_designer!(microgrid, bar(), ω_d)

# Simulate
simulate!(microgrid, controller, designer, ω_a)
