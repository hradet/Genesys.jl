using DataFrames, Genesys, Distributions, CSV, JLD, Dates

# Load data from CSV
data_op = CSV.read(joinpath(pwd(),"data","data_op.csv"), dateformat = "d/m/y HH:MM")
data_inv = CSV.read(joinpath(pwd(),"data","data_inv.csv"))

# Initialize distributions parameters - If the mean value is requested, set the variance to zero...
parameters = (
Δh = 1, # in hours
H = 365 * 24, # in hours
Δy = 1, # in years
Y = 20, # in years
ns = 2 * 1, # number of scenarios - 1 scenario corresponds to Y years - The first ns/2 scenario will be for optimization and the rest for simulation
)

scenarios = (
timestamp = data_op.timestamp,
ld_E = (n_markov = 7, data = data_op.ld_E),
ld_H = (n_markov = 7, data = data_op.ld_H),
pv_E = (n_markov = 10, data = data_op.pv),
C_pv = (distribution = Normal.(data_inv.C_pv, 0.07 * data_inv.C_pv), data = data_inv.C_pv),
C_liion = (distribution = Normal.(data_inv.C_liion, 0.07 * data_inv.C_liion), data = data_inv.C_liion),
C_tes = (distribution = Normal.(data_inv.C_tes, 0.07 * data_inv.C_tes), data = data_inv.C_tes),
C_tank = (distribution = Normal.(data_inv.C_tank, 0.07 * data_inv.C_tank), data = data_inv.C_tank),
C_elyz = (distribution = Normal.(data_inv.C_elyz, 0.07 * data_inv.C_elyz), data = data_inv.C_elyz),
C_fc = (distribution = Normal.(data_inv.C_fc, 0.07 * data_inv.C_fc), data = data_inv.C_fc),
C_heater = (distribution = Normal.(data_inv.C_heater, 0.07 * data_inv.C_heater), data = data_inv.C_heater),
C_grid_out = data_op.C_grid_out,
C_grid_in = data_op.C_grid_in,
)

# Store the parameters in a dict
outputGUI = Dict(
"scenarios" => scenarios,
"parameters" => parameters
)

# Initialize generation and demand profiles scenarios
pv_E, ld_E, ld_H = Genesys.initialize_generation_demand_scenario(outputGUI)

# Initialize investment cost scenarios
C_pv, C_liion, C_tes, C_tank, C_elyz, C_fc, C_heater = Genesys.initialize_investment_cost_scenario(outputGUI)

# Inititalize operating cost scenarios
C_grid_in, C_grid_out = Genesys.initialize_operating_cost_scenario(outputGUI)

# Even number scenarios for optimization...
ω_optim = Dict(
"timestamp" => Vector(scenarios.timestamp),
"pv_E" => pv_E[:,:,2:2:end],
"ld_E" => ld_E[:,:,2:2:end],
"ld_H" => ld_H[:,:,2:2:end],
"C_pv" => C_pv[:,2:2:end],
"C_liion" => C_liion[:,2:2:end],
"C_tes" => C_tes[:,2:2:end],
"C_tank" => C_tank[:,2:2:end],
"C_elyz" => C_elyz[:,2:2:end],
"C_fc" => C_fc[:,2:2:end],
"C_heater" => C_heater[:,2:2:end],
"C_grid_in" => C_grid_in[:,:,2:2:end],
"C_grid_out" => C_grid_out[:,:,2:2:end],
)

# ... and odd number scenarios for simulation
ω_simu = Dict(
"timestamp" => Vector(scenarios.timestamp),
"pv_E" => pv_E[:,:,1:2:end],
"ld_E" => ld_E[:,:,1:2:end],
"ld_H" => ld_H[:,:,1:2:end],
"C_pv" => C_pv[:,1:2:end],
"C_liion" => C_liion[:,1:2:end],
"C_tes" => C_tes[:,1:2:end],
"C_tank" => C_tank[:,1:2:end],
"C_elyz" => C_elyz[:,1:2:end],
"C_fc" => C_fc[:,1:2:end],
"C_heater" => C_heater[:,1:2:end],
"C_grid_in" => C_grid_in[:,:,1:2:end],
"C_grid_out" => C_grid_out[:,:,1:2:end],
)

# Save scenarios as .jld files
save(joinpath("data", "input_data_stochastic.jld"), "ω_optim", ω_optim, "ω_simu", ω_simu)

#= Command to load the .jld file :
 ω_optim, ω_simu = load(joinpath("data", "input_data.jld"), "ω_optim", "ω_simu")
 =#
