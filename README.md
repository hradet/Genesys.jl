# Genesys

A generic module written in Julia to asses and compare different design and control approaches for distributed energy systems (DES) in a stochastic framework. The `simulate!` function includes multi-stage investment periods and multi-scenarios assessment.  

Note that **_out-of-sample_** assesment is made possible as we can optimize and simulate with different scenarios.

# Installation
In order to use the package, follow the [managing package guideline](https://julialang.github.io/Pkg.jl/v1/managing-packages/) for uneregistred packages.

# Resolutions methods
- Design
  - Manual
  - MILP 
  - Metaheuristic
 
- Control
  - Dummy
  - Anticipative
  - Rule based (RBC)
  - Open Loop Feedback Control (OLFC) - OLFC with a single scenario is equivalent to MPC...
  
# Example
We provide a simple example with the dummy controller and designer.

```Julia
# Load packages
using Genesys, CSV, DataFrames, JLD, Dates

# Constant
const nh, ny, ns = 8760, 20, 1 # nh = operation stages, ny = investment stages, ns = scenarios

# Load data
data = load(joinpath("data","ausgrid_scenarios.jld"))

# Initialize scenarios
ω_optim, ω_simu = Scenarios(data["ω_optim"], 1:nh, 1:ny, 1:ns), Scenarios(data["ω_simu"],  1:nh, 1:ny, 1:ns)

# Initialize DES
DES = DistributedEnergySystem(ld_E = Load(),
                              pv = Source(),
                              liion = Liion(),                             
                              grid = Grid(),
                              parameters = Genesys.GlobalParameters(nh, ny, ns, τ_share = 0.8))

# Initialize controller
controller = initialize_controller!(DES, Dummy(), ω_optim)

# Initialize designer
designer = initialize_designer!(DES, Manual(), ω_optim)

# Simulate
simulate!(DES, controller, designer, ω_simu, options = Genesys.Options(mode="multithreads"))

```
