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
mutable struct DummyController <: Genesys.AbstractController
    u
    DummyController() = new()
end
mutable struct DummyDesigner <: Genesys.AbstractDesigner
    u
    DummyDesigner() = new()
end

# Define offline computations
function Genesys.initialize_designer(ld::Genesys.Load, pv::Genesys.Source,
    liion::Genesys.Liion, designer::DummyDesigner,
     grid::Genesys.Grid, ω_optim::Genesys.Scenarios, parameters::NamedTuple)
     # Parameters
     nh = size(ld.power_E,1) # number of simulation hours in one year
     ny = size(ld.power_E,2) # number of simulation years
     ns = size(ld.power_E,3) # number of scenarios

     # Initialize controller and designer policies
     # The policy must be initialize at this place...

     # Initialize decisions variables
     designer.u = (
     u_liion = zeros(ny,ns),
     u_pv = zeros(ny,ns),
     )
end

function Genesys.initialize_controller(ld::Genesys.Load, pv::Genesys.Source,
    liion::Genesys.Liion, controller::DummyController,
     grid::Genesys.Grid, ω_optim::Genesys.Scenarios, parameters::NamedTuple)
     # Parameters
     nh = size(ld.power_E,1) # number of simulation hours in one year
     ny = size(ld.power_E,2) # number of simulation years
     ns = size(ld.power_E,3) # number of scenarios

     # Initialize controller and designer policies
     # The policy must be initialize at this place...

     # Initialize decisions variables
     controller.u = (
     u_liion = zeros(nh,ny,ns),
     )
end

# Define online functions to compute decisions
function Genesys.compute_operation_decisions(h::Int64, y::Int64, s::Int64, ld::Genesys.Load,
    pv::Genesys.Source, liion::Genesys.Liion, grid::Genesys.Grid, controller::DummyController,
    ω_optim::Genesys.Scenarios, parameters::NamedTuple)
    return 0
end
function Genesys.compute_investment_decisions(y::Int64, s::Int64, ld::Genesys.Load, pv::Genesys.Source,
     liion::Genesys.Liion ,grid::Genesys.Grid,
     designer::DummyDesigner, ω_optim::Genesys.Scenarios, parameters::NamedTuple)
    return 0, 0
end

#=
Let's simulate the microgrid with the dummies controller and designer...
=#

# GUI loading
outputGUI = Genesys.loadGUI("", "")

# Initialization without any controller and designer
ld, pv, liion, _, _, _, _, _, _, _, grid, ω_optim, ω_simu = Genesys.initialization(outputGUI)

# Intialize the controller and designer
dummy_controller, dummy_designer = DummyController(), DummyDesigner()

# Initialize designer
Genesys.initialize_designer(ld, pv, liion, dummy_designer, grid, ω_optim, outputGUI["parameters"])

# Initialize controller
Genesys.initialize_controller(ld, pv, liion, dummy_controller, grid, ω_optim, outputGUI["parameters"])

# Simulate
Genesys.simulate(ld, pv, liion, dummy_controller, dummy_designer, grid, ω_optim, ω_simu, outputGUI["parameters"])
