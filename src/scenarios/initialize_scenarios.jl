# Initialize generation and demand scenarios
function initialize_generation_demand_scenario(outputGUI)
    # Parameters
    nh = size(outputGUI["scenarios"].ld_E.data,1) # number of step for one year scenario
    ny = length(outputGUI["parameters"].Δy:outputGUI["parameters"].Δy:outputGUI["parameters"].Y)
    ns = outputGUI["parameters"].ns # number of ny years scenarios

    # Pre-allocate
    pv_E = zeros(nh, ns * ny)
    ld_E = zeros(nh, ns * ny)
    ld_H = zeros(nh, ns * ny)

    # Clustering data
    # PV
    pv_E_month = Genesys.clustering_month(outputGUI["scenarios"].pv_E.data, outputGUI["scenarios"].timestamp)
    # Elec
    ld_E_wk = Genesys.clustering_month_week(outputGUI["scenarios"].ld_E.data, outputGUI["scenarios"].timestamp)
    ld_E_wkd = Genesys.clustering_month_weekend(outputGUI["scenarios"].ld_E.data, outputGUI["scenarios"].timestamp)
    # Heat
    ld_H_wk = Genesys.clustering_month_week(outputGUI["scenarios"].ld_H.data, outputGUI["scenarios"].timestamp)
    ld_H_wkd = Genesys.clustering_month_weekend(outputGUI["scenarios"].ld_H.data, outputGUI["scenarios"].timestamp)

    # Compute markov chain
    # PV
    mc_pv_E = Genesys.compute_markovchain(pv_E_month, outputGUI["scenarios"].pv_E.n_markov)
    # Elec
    mc_ld_E_wk = Genesys.compute_markovchain(ld_E_wk, outputGUI["scenarios"].ld_E.n_markov)
    mc_ld_E_wkd = Genesys.compute_markovchain(ld_E_wkd, outputGUI["scenarios"].ld_E.n_markov)
    # Heat
    mc_ld_H_wk = Genesys.compute_markovchain(ld_H_wk, outputGUI["scenarios"].ld_H.n_markov)
    mc_ld_H_wkd = Genesys.compute_markovchain(ld_H_wkd, outputGUI["scenarios"].ld_H.n_markov)

    # Compute scenarios from MC
    @showprogress for s in 1:ns * ny
        pv_E[:,s] = Genesys.compute_scenario(mc_pv_E, outputGUI["scenarios"].pv_E.data[1], outputGUI["scenarios"].timestamp[1], 1, nh-1)
        ld_E[:,s] = Genesys.compute_scenario(mc_ld_E_wk, mc_ld_E_wkd, outputGUI["scenarios"].ld_E.data[1], outputGUI["scenarios"].timestamp[1], 1, nh-1)
        ld_H[:,s] = Genesys.compute_scenario(mc_ld_H_wk, mc_ld_H_wkd, outputGUI["scenarios"].ld_H.data[1], outputGUI["scenarios"].timestamp[1], 1, nh-1)
    end

    # Reshape to ny years scenarios
    pv_E = reshape(pv_E, nh, ny, ns)
    ld_E = reshape(ld_E, nh, ny, ns)
    ld_H = reshape(ld_H, nh, ny, ns)

    return pv_E, ld_E, ld_H
end
# Inititalize investment cost scenarios
function initialize_investment_cost_scenario(outputGUI)
    # Parameters
    ns = outputGUI["parameters"].ns # number of scenarios
    # Initialize investment cost scenarios
    C_pv = compute_scenarios(outputGUI["scenarios"].C_pv, ns)
    C_liion = compute_scenarios(outputGUI["scenarios"].C_liion, ns)
    C_tes = compute_scenarios(outputGUI["scenarios"].C_tes, ns)
    C_tank = compute_scenarios(outputGUI["scenarios"].C_tank, ns)
    C_elyz = compute_scenarios(outputGUI["scenarios"].C_elyz, ns)
    C_fc = compute_scenarios(outputGUI["scenarios"].C_fc, ns)
    C_heater = compute_scenarios(outputGUI["scenarios"].C_heater, ns)

    return C_pv, C_liion, C_tes, C_tank, C_elyz, C_fc, C_heater
end
# Compute scenario from distribution
function compute_scenarios(scenario, ns)
    # Pre-allocate
    scenarios = Array{Float64,2}(undef, size(scenario.data,1), ns)
    # Generate scenarios
    for (y, distribution) in enumerate(scenario.distribution)
        scenarios[y,:] = rand(distribution, ns)
    end

    return scenarios
end
# Inititalize operating cost scenarios
function initialize_operating_cost_scenario(outputGUI)
    # Parameters
    nh = size(outputGUI["scenarios"].ld_E.data,1) # number of step for one year scenario
    ny = length(outputGUI["parameters"].Δy:outputGUI["parameters"].Δy:outputGUI["parameters"].Y)
    ns = outputGUI["parameters"].ns # number of ny years scenarios

    # Repeat the values ns * ny times
    C_grid_in = repeat(outputGUI["scenarios"].C_grid_in, 1, ny, ns)
    C_grid_out = repeat(outputGUI["scenarios"].C_grid_out, 1, ny, ns)

    return C_grid_in, C_grid_out
end
