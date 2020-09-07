
# initialization
function initialization(outputGUI)
    #---------------------------------------------------------------------------
    #                                  Parameters
    #---------------------------------------------------------------------------
    nh = min(length(outputGUI["parameters"].Δh:outputGUI["parameters"].Δh:outputGUI["parameters"].H), size(outputGUI["scenarios"].ω_optim["ld_E"],1))
    ny = min(length(outputGUI["parameters"].Δy:outputGUI["parameters"].Δy:outputGUI["parameters"].Y), size(outputGUI["scenarios"].ω_optim["ld_E"],2))
    ns = min(outputGUI["parameters"].ns, size(outputGUI["scenarios"].ω_optim["ld_E"],3))

    # TODO: Modifier ns et Scenarios !
    #---------------------------------------------------------------------------
    #                                  Scenarios
    #---------------------------------------------------------------------------
    ω_optim = Scenarios(outputGUI["scenarios"].ω_optim, nh, ny, ns)
    ω_simu = Scenarios(outputGUI["scenarios"].ω_simu, nh, ny, ns)

    #---------------------------------------------------------------------------
    #                                  Load
    #---------------------------------------------------------------------------
    ld=Load(outputGUI["ld"], nh, ny, ns)

    #---------------------------------------------------------------------------
    #                                  Source
    #---------------------------------------------------------------------------
    pv=Source(outputGUI["pv"], nh, ny, ns)
    pv.powerMax[1,:] .= outputGUI["pv"].powerMax_ini

    #---------------------------------------------------------------------------
    #                                  Li-ion
    #---------------------------------------------------------------------------
    liion = Liion(outputGUI["liion"], nh, ny, ns)
    liion.Erated[1,:] .= outputGUI["liion"].Erated_ini
    liion.soc[1,1,:] .= outputGUI["liion"].α_soc_ini
    liion.soh[1,1,:] .= outputGUI["liion"].soh_ini

    #---------------------------------------------------------------------------
    #                                  Electrolyzer
    #---------------------------------------------------------------------------
    elyz = Electrolyzer(outputGUI["elyz"], nh, ny, ns)
    elyz.powerMax[1,:] .= outputGUI["elyz"].powerMax_ini
    elyz.soh[1,1,:] .= outputGUI["elyz"].soh_ini

    #---------------------------------------------------------------------------
    #                                  Fuel Cell
    #---------------------------------------------------------------------------
    fc = FuelCell(outputGUI["fc"], nh, ny, ns)
    fc.powerMax[1,:] .= outputGUI["fc"].powerMax_ini
    fc.soh[1,1,:] .= outputGUI["fc"].soh_ini

    #---------------------------------------------------------------------------
    #                                  H2 Tank
    #---------------------------------------------------------------------------
    h2tank = H2Tank(outputGUI["h2tank"], nh, ny, ns)
    h2tank.Erated[1,:] .= outputGUI["h2tank"].Erated_ini
    h2tank.soc[1,1,:] .= outputGUI["h2tank"].α_soc_ini

    #---------------------------------------------------------------------------
    #                                  Thermal Energy Storage
    #---------------------------------------------------------------------------
    tes = ThermalSto(outputGUI["tes"], nh, ny, ns)
    tes.Erated[1,:] .= outputGUI["tes"].Erated_ini
    tes.soc[1,1,:] .= outputGUI["tes"].α_soc_ini

    #---------------------------------------------------------------------------
    #                                  Heater
    #---------------------------------------------------------------------------
    heater = Heater(outputGUI["heater"], nh, ny, ns)
    heater.powerMax .= outputGUI["heater"].powerMax_ini

    #---------------------------------------------------------------------------
    #                                  Grid
    #---------------------------------------------------------------------------
    grid = Grid(outputGUI["grid"], nh, ny, ns)

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
