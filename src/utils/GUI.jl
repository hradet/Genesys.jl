# Load parameters from virtual GUI
function loadGUI(controller_flag::String, designer_flag::String; ns = 1)
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
    τ_share = 0., # share of renwables
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
    soc_ini = 0.5,
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
    soc_ini = 0.5,
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
    soc_ini = 0.5,
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
    eac_stoch = Dict(
        "reduction" => "manual",
        "idx_scenario" => 1,
        "idx_years" => 1:20,
        "risk" => "esperance"),
    eac_stoch_td = Dict(),
    metaheuristic = Dict(
        "method" => "cmaes",
        "reduction" => "manual",
        "idx_scenario" => 1,
        "idx_year" => 1,
        "u0" => (liion = 0., pv = 0.),
        "lb" => (liion = 0., pv = 0.),
        "ub" => (liion = 1000., pv = 1000.),
        "constraints" => (λ = 1e10,),
        "iteration" => 10,
        "reltol" => 0.05),
    )

    #---------------------------------------------------------------------------
    #                                  Grid
    #---------------------------------------------------------------------------
    grid = (
    powerMax = 36.,
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

# Initialization from GUI
function initialization(outputGUI)
    #---------------------------------------------------------------------------
    #                                  Parameters
    #---------------------------------------------------------------------------
    nh = min(length(outputGUI["parameters"].Δh:outputGUI["parameters"].Δh:outputGUI["parameters"].H), size(outputGUI["scenarios"].ω_optim["ld_E"],1))
    ny = min(length(outputGUI["parameters"].Δy:outputGUI["parameters"].Δy:outputGUI["parameters"].Y), size(outputGUI["scenarios"].ω_optim["ld_E"],2))
    ns = min(outputGUI["parameters"].ns, size(outputGUI["scenarios"].ω_optim["ld_E"],3))

    #---------------------------------------------------------------------------
    #                                  Scenarios
    #---------------------------------------------------------------------------
    ω_optim = Scenarios(outputGUI["scenarios"].ω_optim, nh, ny, ns)
    ω_simu = Scenarios(outputGUI["scenarios"].ω_simu, nh, ny, ns)

    #---------------------------------------------------------------------------
    #                                  Load
    #---------------------------------------------------------------------------
    ld = Load()

    # Preallocation
    preallocate!(ld, nh, ny, ns)

    #---------------------------------------------------------------------------
    #                                  Source
    #---------------------------------------------------------------------------
    pv = Source(lifetime = outputGUI["pv"].lifetime,
                powerMax_ini = outputGUI["pv"].powerMax_ini)

    # Preallocation
    preallocate!(pv, nh, ny, ns)

    #---------------------------------------------------------------------------
    #                                  Li-ion
    #---------------------------------------------------------------------------
    liion = Liion(α_p_ch = outputGUI["liion"].α_p_ch,
            α_p_dch = outputGUI["liion"].α_p_dch,
            η_ch = outputGUI["liion"].η_ch,
            η_dch = outputGUI["liion"].η_dch,
            η_self = outputGUI["liion"].η_self,
            α_soc_min = outputGUI["liion"].α_soc_min,
            α_soc_max = outputGUI["liion"].α_soc_max,
            lifetime = outputGUI["liion"].lifetime,
            nCycle = outputGUI["liion"].nCycle,
            dod = outputGUI["liion"].dod,
            Erated_ini = outputGUI["liion"].Erated_ini,
            soh_ini = outputGUI["liion"].soh_ini,
            soc_ini = outputGUI["liion"].soc_ini)

    # Preallocation
    preallocate!(liion, nh, ny, ns)

    #---------------------------------------------------------------------------
    #                                  Electrolyzer
    #---------------------------------------------------------------------------
    elyz = Electrolyzer(α_p = outputGUI["elyz"].α_p,
                η_E_H2 = outputGUI["elyz"].η_E_H2,
                η_E_H = outputGUI["elyz"].η_E_H,
                lifetime = outputGUI["elyz"].lifetime,
                nHoursMax = outputGUI["elyz"].nHoursMax,
                powerMax_ini = outputGUI["elyz"].powerMax_ini,
                soh_ini = outputGUI["elyz"].soh_ini)

    # Preallocation
    preallocate!(elyz, nh, ny, ns)

    #---------------------------------------------------------------------------
    #                                  Fuel Cell
    #---------------------------------------------------------------------------
    fc = FuelCell(α_p = outputGUI["fc"].α_p,
                η_H2_E = outputGUI["fc"].η_H2_E,
                η_H2_H = outputGUI["fc"].η_H2_H,
                lifetime = outputGUI["fc"].lifetime,
                nHoursMax = outputGUI["fc"].nHoursMax,
                powerMax_ini = outputGUI["fc"].powerMax_ini,
                soh_ini = outputGUI["fc"].soh_ini)

    # Preallocation
    preallocate!(fc, nh, ny, ns)

    #---------------------------------------------------------------------------
    #                                  H2 Tank
    #---------------------------------------------------------------------------
    h2tank = H2Tank(α_p_ch = outputGUI["h2tank"].α_p_ch,
            α_p_dch = outputGUI["h2tank"].α_p_dch,
            η_ch = outputGUI["h2tank"].η_ch,
            η_dch = outputGUI["h2tank"].η_dch,
            η_self = outputGUI["h2tank"].η_self,
            α_soc_min = outputGUI["h2tank"].α_soc_min,
            α_soc_max = outputGUI["h2tank"].α_soc_max,
            lifetime = outputGUI["h2tank"].lifetime,
            Erated_ini = outputGUI["h2tank"].Erated_ini,
            soc_ini = outputGUI["h2tank"].soc_ini)

    # Preallocation
    preallocate!(h2tank, nh, ny, ns)

    #---------------------------------------------------------------------------
    #                                  Thermal Energy Storage
    #---------------------------------------------------------------------------
    tes = ThermalSto(α_p_ch = outputGUI["tes"].α_p_ch,
            α_p_dch = outputGUI["tes"].α_p_dch,
            η_ch = outputGUI["tes"].η_ch,
            η_dch = outputGUI["tes"].η_dch,
            η_self = outputGUI["tes"].η_self,
            α_soc_min = outputGUI["tes"].α_soc_min,
            α_soc_max = outputGUI["tes"].α_soc_max,
            lifetime = outputGUI["tes"].lifetime,
            Erated_ini = outputGUI["tes"].Erated_ini,
            soc_ini = outputGUI["tes"].soc_ini)

    # Preallocation
    preallocate!(tes, nh, ny, ns)

    #---------------------------------------------------------------------------
    #                                  Heater
    #---------------------------------------------------------------------------
    heater = Heater(η_E_H = outputGUI["heater"].η_E_H,
                    powerMax_ini = outputGUI["heater"].powerMax_ini,
                    lifetime = outputGUI["heater"].lifetime)

    # Preallocation
    preallocate!(heater, nh, ny, ns)

    #---------------------------------------------------------------------------
    #                                  Grid
    #---------------------------------------------------------------------------
    grid = Grid(powerMax = outputGUI["grid"].powerMax)

    # Preallocation
    preallocate!(grid, nh , ny, ns)

    #---------------------------------------------------------------------------
    #                                  Controller
    #---------------------------------------------------------------------------
    if outputGUI["controller"].id == "AnticipativeController"
        controller = AnticipativeController()
    elseif outputGUI["controller"].id == "RuleBasedController"
        controller = RuleBasedController()
        controller.parameters = outputGUI["controller"].rb
    elseif outputGUI["controller"].id == "MPCController"
        controller = MPCController()
        controller.parameters = outputGUI["controller"].mpc
    else
        controller = nothing
        println("Unknown controller id...")
    end

    #---------------------------------------------------------------------------
    #                                  Designer
    #---------------------------------------------------------------------------
    if outputGUI["designer"].id == "AnticipativeMultiStageDesigner"
        designer = AnticipativeMultiStageDesigner()
        designer.parameters = outputGUI["designer"].anticipative_multi
    elseif outputGUI["designer"].id == "AnticipativeOneStageDesigner"
        designer = AnticipativeOneStageDesigner()
        designer.parameters = outputGUI["designer"].anticipative_one
    elseif outputGUI["designer"].id == "AnticipativeOneStageOnlineDesigner"
        designer = AnticipativeOneStageOnlineDesigner()
        designer.parameters = outputGUI["designer"].anticipative_one_up
    elseif outputGUI["designer"].id == "DummyDesigner"
        designer = DummyDesigner()
    elseif outputGUI["designer"].id == "RuleBasedDesigner"
        designer = RuleBasedDesigner()
        designer.parameters = outputGUI["designer"].rb
    elseif outputGUI["designer"].id == "EACDesigner"
        designer = EACDesigner()
        designer.parameters = outputGUI["designer"].eac
    elseif outputGUI["designer"].id == "EACStochasticDesigner"
        designer = EACStochasticDesigner()
        designer.parameters = outputGUI["designer"].eac_stoch
    elseif outputGUI["designer"].id == "TypicalDayEACDesigner"
        designer = TypicalDayEACDesigner()
        designer.parameters = outputGUI["designer"].eac_td
    elseif outputGUI["designer"].id == "MetaHeuristicDesigner"
        designer = MetaHeuristicDesigner()
        designer.parameters = outputGUI["designer"].metaheuristic
    else
        designer = nothing
        println("Unknown designer id...")
    end

   return ld, pv, liion, h2tank, elyz, fc, tes, heater, controller, designer, grid, ω_optim, ω_simu
end
