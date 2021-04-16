# Genesys

A generic module written in Julia to asses and compare different design and control methods for distributed energy systems (DES) in a stochastic framework. The `simulate!` function includes multi-stage investment periods and multi-scenarios assessment.  

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
  - Stochastic dual dynamic programming (SDDP)
 
 # Scenario generation and reduction methods
- Generation  
  - Anticipative
  - Markov

- Reduction
  - Manual
  - Mean value
  - Monte Carlo  
  - Feature based
  
# Example
We provide a simple example with the rule-based controller and metaheuristic designer.

```Julia
# Load packages
using Genesys, JLD, Dates

# Constant
const nh, ny, ns = 8760, 2, 1000 # nh = operation stages, ny = investment stages, ns = scenarios

# Load data
data = load("...path\scenarios.jld")

# Initialize scenarios
ω_optim, ω_simu = Scenarios(data["ω_optim"], 1:nh, 1:ny, 1:ns), Scenarios(data["ω_simu"],  1:nh, 1:ny, 1:ns)

# Initialize DES
DES = DistributedEnergySystem(ld_E = Load(),
                              pv = Source(),
                              liion = Liion(),                             
                              grid = Grid(),
                              parameters = Genesys.GlobalParameters(nh, ny, ns, renewable_share = 0.8))

# Initialize controller
controller = initialize_controller!(DES, RBC(), ω_optim)

# Initialize designer
designer = initialize_designer!(DES, Metaheuristic(), ω_optim)

# Simulate
simulate!(DES, controller, designer, ω_simu, options = Genesys.Options(mode="multithreads"))

# Metrics
metrics = Metrics(DES, designer)


```
