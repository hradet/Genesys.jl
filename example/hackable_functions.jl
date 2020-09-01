# Genesys
using Genesys
# Optimization
using JuMP
# Lecture CSV
using CSV, DataFrames, JLD2
# Plot
using Seaborn
pygui(true)

# GUI function
include("load_GUI.jl")

#=
You can either use a controller and designer from the library,
or you can define your own controller and designer and modify the associated functions...
=#

# Define your own controller and designer
mutable struct DummyController <: Genesys.AbstractController
    π::Function
    u::NamedTuple
    DummyController() = new()
end
mutable struct DummyDesigner <: Genesys.AbstractDesigner
    π::Function
    u::NamedTuple
    DummyDesigner() = new()
end

# Define offline computations
function Genesys.offline_optimization(ld::Genesys.Load, pv::Genesys.Source,
    liion::Genesys.Liion, controller::DummyController, designer::DummyDesigner,
     grid::Genesys.Grid, ω_optim::Genesys.Scenarios, parameters::NamedTuple)
     # Parameters
     nh = size(ld.power_E,1) # number of simulation hours in one year
     ny = size(ld.power_E,2) # number of simulation years
     ns = size(ld.power_E,3) # number of scenarios

     # Initialize controller and designer policies
     # The policy must be initialize at this place...

     # Initialize decisions variables
     # Operation decisions
     controller.u = (
     u_liion = zeros(nh,ny,ns),
     )
     # Investment controls
     designer.u = (
     u_liion = zeros(ny,ns),
     u_pv = zeros(ny,ns),
     )
end

# Define online functions to compute decisions
function Genesys.compute_operation_decisions(h::Int64, y::Int64, s::Int64, ld::Genesys.Load,
    pv::Genesys.Source, liion::Genesys.Liion, grid::Genesys.Grid, controller::DummyController,
    ω_optim::Genesys.Scenarios, parameters::NamedTuple)
    return 0
end
function Genesys.compute_investment_decisions(y::Int64, s::Int64, ld::Genesys.Load, pv::Genesys.Source,
     liion::Genesys.Liion ,grid::Genesys.Grid, controller::DummyController,
     designer::DummyDesigner, ω_optim::Genesys.Scenarios, parameters::NamedTuple)
    return 0, 0
end

#=
Let's simulate the microgrid with the dummies controller and designer...
=#

# GUI loading
outputGUI = loadGUI("", "", 1., 0., 0.)

# Initialization without any controller and designer
ld, pv, liion, _, _, _, _, _, _, _, grid, ω_optim, ω_simu = Genesys.initialization(outputGUI)

# Intialize the controller and designer
dummy_controller, dummy_designer = DummyController(), DummyDesigner()

# Offline phase
Genesys.offline_optimization(ld, pv, liion, dummy_controller, dummy_designer, grid, ω_optim, outputGUI["parameters"])

# Simulate
costs = Genesys.simulate(ld, pv, liion, dummy_controller, dummy_designer, grid, ω_optim, ω_simu, outputGUI["parameters"])
