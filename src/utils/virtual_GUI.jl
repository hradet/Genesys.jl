# Load parameters from GUI
function loadGUI(controller_flag::String, designer_flag::String; τ_power = 0., τ_energy = 0., ns = 1)
    # Directory path
    dirPath=pwd()

    #---------------------------------------------------------------------------
    #                          Global parameters
    #---------------------------------------------------------------------------
    parameters = (
    Δh = 1, # in hours
    H = 365*24, # in hours
    Δy = 1, # in years
    Y = 20, # in years
    τ = 0.045, # discount rate
    ns = ns, # number of scenarios of ny years
    )

    #---------------------------------------------------------------------------
    #                                  loading data from JLD
    #---------------------------------------------------------------------------
    # TODO filename en input de la fonction
    scenarios = (
    ω_optim = load(joinpath("data", "input_data_stochastic.jld"), "ω_optim"),
    ω_simu = load(joinpath("data", "input_data_stochastic.jld"), "ω_simu"),
    )

    #---------------------------------------------------------------------------
    #                                  Load
    #---------------------------------------------------------------------------
    ld = ()

    #---------------------------------------------------------------------------
    #                                  Source
    #---------------------------------------------------------------------------
    pv = (
    lifetime = 25.,
    powerMax_ini = 0.,
    )

    #---------------------------------------------------------------------------
    #                                  Li-ion
    #---------------------------------------------------------------------------
    liion = (
    Erated_ini = 0.,
    α_p_ch = 1.5,
    α_p_dch = 1.5,
    η_ch = 0.8,
    η_dch = 0.8,
    η_self = 0.,
    α_soc_min = 0.2,
    α_soc_max = 0.8,
    α_soc_ini = 0.5,
    # Aging
    lifetime = 12.,
    soh_ini=1.,
    nCycle = 2500,#7000,
    dod = 0.6,
    )

    #---------------------------------------------------------------------------
    #                                  Electrolyzer
    #---------------------------------------------------------------------------
    elyz = (
    powerMax_ini = 0.,
    α_p = 5/100,
    η_E_H2 = 0.5,
    η_E_H = 0.3,
    # Aging
    lifetime = 10.,
    soh_ini = 1.,
    nHoursMax = 26000,
    )

    #---------------------------------------------------------------------------
    #                                  Fuel Cell
    #---------------------------------------------------------------------------
    fc = (
    powerMax_ini = 0.,
    α_p = 8/100,
    η_H2_E = 0.4,
    η_H2_H = 0.4,
    # Aging
    lifetime = 5.,
    soh_ini = 1.,
    nHoursMax = 10000,
    )

    #---------------------------------------------------------------------------
    #                                  H2 Tank
    #---------------------------------------------------------------------------
    h2tank =(
    Erated_ini = 0.,
    α_p_ch = 1000.,
    α_p_dch = 1000.,
    η_ch = 1.,
    η_dch = 1.,
    η_self = 0.,
    α_soc_min = 0.,
    α_soc_max = 1.,
    α_soc_ini = 0.5,
    lifetime = 25.,
    )

    #---------------------------------------------------------------------------
    #                                  Thermal Energy Storage
    #---------------------------------------------------------------------------
    tes = (
    Erated_ini = 0.,
    α_p_ch = 1.5,
    α_p_dch = 1.5,
    η_ch = 0.9,
    η_dch = 0.9,
    η_self = 0.01,
    α_soc_min = 0.,
    α_soc_max = 1.,
    α_soc_ini = 0.5,
    lifetime = 25.,
    )

    #---------------------------------------------------------------------------
    #                                  Heater
    #---------------------------------------------------------------------------
    heater = (
    powerMax_ini = maximum(scenarios.ω_optim["ld_H"]) + 5,
    η_E_H = 1.,
    lifetime = 25.,
    )

    #---------------------------------------------------------------------------
    #                                  Controller
    #---------------------------------------------------------------------------
    controller=(
    id = controller_flag,
    mpc = Dict(
        "horizon" => 12), # hours
    rb = Dict(
        "β_min_tes" => 0.2,
        "β_max_tes" => 0.9,
        "β_min_fc" => 0.25,
        "β_max_fc" => 0.3,
        "β_min_elyz" => 0.4,
        "β_max_elyz" => 0.45),

    )

    #---------------------------------------------------------------------------
    #                                  Designer
    #---------------------------------------------------------------------------
    designer=(
    id = designer_flag,
    anticipative_multi = Dict(
        "reduction" => "manual",
        "idx_scenario" => 1),
    anticipative_one = Dict(
        "reduction" => "manual",
        "idx_scenario" => 1,
        "idx_years" => 1:20,
        "risk" => "esperance"),
    anticipative_one_up = Dict(
        "reduction" => "manual",
        "idx_scenario" => 1,
        "idx_years" => 1:20,
        "risk" => "esperance"),
    eac = Dict(
        "reduction" => "manual",
        "idx_scenario" => 1,
        "idx_year" => 1),
    eac_td = Dict(
        "reduction" => "manual",
        "idx_scenario" => 1,
        "idx_year" => 1,
        "ntd" => 5),
    eac_stoch = Dict(
        "reduction" => "manual",
        "idx_scenario" => 1,
        "idx_years" => 1:20,
        "risk" => "esperance"),
    eac_stoch_td = Dict(),
    metaheuristic = Dict(),
    metaheuristic_stoch = Dict(),
    rb = Dict(),
    )

    #---------------------------------------------------------------------------
    #                                  Grid
    #---------------------------------------------------------------------------
    grid = (
    τ_power = τ_power,
    τ_energy = τ_energy,
    )

    #---------------------------------------------------------------------------
    #                                  Output GUI
    #---------------------------------------------------------------------------
    outputGUI = Dict(
    "ld" => ld,
    "pv" => pv,
    "liion" => liion,
    "elyz" => elyz,
    "fc" => fc,
    "h2tank" => h2tank,
    "tes" => tes,
    "heater" => heater,
    "controller" => controller,
    "designer" => designer,
    "grid" => grid,
    "scenarios" => scenarios,
    "parameters" => parameters
    )

    return outputGUI
end
