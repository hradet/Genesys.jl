# Genesys

A generic module written in Julia to asses and compare different design and control approaches for distributed energy systems (DES).

# Installation
In order to use the package, follow the [managing package guideline](https://julialang.github.io/Pkg.jl/v1/managing-packages/) for uneregistred packages.

# Algorithms
- Design
  - Dummy
  - Equivalent annual cost (EAC)
  - Stochastic equivalent annual cost (EACStoch)
  - Metaheuristic
 
- Control
  - Dummy
  - Rule based (RBC)
  - Model predictive control (MPC)
  
# Example
We provide a simple example with the dummy controller and designer.

```Julia
# Load packages
using Genesys, CSV, DataFrames, JLD, Dates

# Constant
const nh, ny, ns = 8760, 20, 1

# Load data
data = load(joinpath("data","input_data_stochastic.jld"))

# Initialize scenarios
ω_optim, ω_simu = Scenarios(data["ω_optim"], nh, ny, ns), Scenarios(data["ω_simu"],  nh, ny, ns)

# Initialize DES
DES = DistributedEnergySystem(ld_E = Load(),
                              ld_H = Load(),
                              pv = Source(),
                              liion = Liion(),
                              tes = ThermalSto(),
                              h2tank = H2Tank(),
                              elyz = Electrolyzer(),
                              fc = FuelCell(),
                              heater = Heater(),
                              grid = Grid(),
                              parameters = Genesys.GlobalParameters(nh = nh, ny = ny, ns = ns, τ_share = 0.8))

# Initialize controller
controller = initialize_controller!(DES, DummyController(), ω_optim)

# Initialize designer
designer = initialize_designer!(DES, DummyDesigner(), ω_optim)

# Simulate
simulate!(DES, controller, designer, ω_simu, options = Genesys.Options(mode="multithreads"))

```
