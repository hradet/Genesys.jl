using CSV, Seaborn, Genesys
pygui(true)

# GUI function
include("load_GUI.jl")

# Parameters
controller_flag, designer_flag = "" , ""
τ_energy_grid, τ_power_grid, τ_cost_grid = 0., 0., 1
ntd, ntb = 50, 2

# Reshape
# GUI loading
outputGUI = loadGUI(controller_flag, designer_flag, τ_cost_grid, τ_power_grid, τ_energy_grid)

# Initialization
_, _, _, _, _, _, _, _, _, _, _, _, ω_simu = Genesys.initialization(outputGUI)

# Clustering time blocks
ω_tb = Genesys.clustering_time_block(ω_simu, ntb)

# Clustering typical days from time blocks
ω_td = Genesys.clustering_typical_day(ω_simu, ntd)

# Build TD profile
ω_build = Genesys.build_td_profile(ω_td, 1, 1, 8760)

sp=subplot(311)
plot(ω_tb.clusters.ld_E[:,1,1], color=(150,150,150)./255)
plot(ω_build.ld_E[:,1,1], color=(150,0,0)./255)
title("Load_E clustering")
subplot(312, sharex=sp)
plot(ω_tb.clusters.ld_H[:,1,1], color=(150,150,150)./255)
plot(ω_build.ld_H[:,1,1], color=(150,0,0)./255)
title("Load_H clustering")
subplot(313, sharex=sp)
plot(ω_tb.clusters.pv_E[:,1,1], color=(150,150,150)./255)
plot(ω_build.pv_E[:,1,1], color=(150,0,0)./255)
title("PV clustering")
