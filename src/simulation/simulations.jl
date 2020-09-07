#=
    This file includes all the funtions needed to simulate the DES
 =#

# Simulate function
# Simple
function simulate(ld::Load, pv::Source, liion::Liion, controller::AbstractController,
     designer::AbstractDesigner, grid::Grid, ω_optim::AbstractScenarios, ω_simu::AbstractScenarios, parameters::NamedTuple;
     mode = "serial")

    # Parameters
    ns = size(ω_simu.values.ld_E,3)

    if mode == "serial"
        # We simulate over the horizon for all the scenarios
        @showprogress for s in 1:ns
            simulate_scenario(s, ld, pv, liion, controller,
                 designer, grid, ω_optim, ω_simu, parameters)
        end
    elseif mode == "multicores"
        # Init
        s = 0
        # We simulate over the horizon for all the scenarios in parallel
        @sync begin
            for p in workers()
                @async begin
                    while true
                        # Increment scenario
                        s += 1
                        # Break if ns is reached...
                        if s > ns
                            break
                        end
                        # Execute the function
                        remotecall_fetch(simulate_scenario, p, s, ld, pv, liion, controller,
                             designer, grid, ω_optim, ω_simu, parameters)
                    end
                end
            end
        end
    elseif mode == "distributed"
        # We simulate over the horizon for all the scenarios in parallel using the distributed macro
        @sync @distributed for s in 1:ns
            simulate_scenario(s, ld, pv, liion, controller,
                 designer, grid, ω_optim, ω_simu, parameters)
        end
    elseif mode == "multithreads"
        # We simulate over the horizon for all the scenarios in parallel using the distributed macro
        Threads.@threads for s in 1:ns
            simulate_scenario(s, ld, pv, liion, controller,
                 designer, grid, ω_optim, ω_simu, parameters)
        end
    else
        println("Unknown mode... Please chose between 'serial', 'multicores', 'distributed' or 'multithreads'.")
    end
end
# Multi-energy
function simulate(ld::Load, pv::Source, liion::Liion, h2tank::H2Tank,
    elyz::Electrolyzer, fc::FuelCell, tes::ThermalSto, heater::Heater,
    controller::AbstractController, designer::AbstractDesigner, grid::Grid,
    ω_optim::AbstractScenarios, ω_simu::AbstractScenarios, parameters::NamedTuple; mode="serial")

    # Parameters
    ns = size(ω_simu.values.ld_E,3)

    # We simulate over the horizon for all the scenarios
    if mode == "serial"
        @showprogress for s in 1:ns
            simulate_scenario(s, ld, pv, liion, h2tank, elyz, fc, tes, heater,
            controller, designer, grid, ω_optim, ω_simu, parameters)
        end
    elseif mode == "multicores"
        #TODO
    elseif mode == "distributed"
        #TODO
    elseif mode == "multithreads"
        #TODO
    else
        println("Unknown mode... Please chose between 'serial', 'multicores', 'distributed' or 'multithreads'.")

    end
end

# Simulate one scenario
# Simple
function simulate_scenario(s::Int64, ld::Load, pv::Source, liion::Liion, controller::AbstractController,
     designer::AbstractDesigner, grid::Grid, ω_optim::AbstractScenarios, ω_simu::AbstractScenarios, parameters::NamedTuple)

    # Parameters
    nh = size(ω_simu.values.ld_E,1)
    ny = size(ω_simu.values.ld_E,2)

    # We simulate over the horizon for a signle scenario
    for y = 1 : ny
            for h = 1 : nh

                # Update operation informations
                update_operation_informations(h, y, s, ld, pv, liion, grid, ω_simu)

                # Compute operation decision variables
                compute_operation_decisions(h, y, s, ld, pv, liion, grid, controller, ω_optim, parameters)

                # Compute operation dynamics
                compute_operation_dynamics(h, y, s, liion, controller, parameters)

                # Power balance constraint checked for each node
                power_balance(h, y, s, ld, pv, liion, grid)

            end

            # Update investment informations
            update_investment_informations(y, s, pv, liion, ω_simu)

            # Compute investment decision variables
            compute_investment_decisions(y, s, ld, pv, liion, grid, designer, ω_optim, parameters)

            # Compute investment dynamics
            compute_investment_dynamics(y, s, pv, liion, designer)

    end
end
# Multi-energy
function simulate_scenario(s::Int64, ld::Load, pv::Source, liion::Liion, h2tank::H2Tank,
    elyz::Electrolyzer, fc::FuelCell, tes::ThermalSto, heater::Heater,
    controller::AbstractController, designer::AbstractDesigner, grid::Grid,
    ω_optim::AbstractScenarios, ω_simu::AbstractScenarios, parameters::NamedTuple)

    # Parameters
    nh = size(ω_simu.values.ld_E,1)
    ny = size(ω_simu.values.ld_E,2)

    # We simulate over the horizon for a signle scenario
    for y = 1 : ny
        for h = 1 : nh

            # Update operation informations
            update_operation_informations(h, y, s, ld, pv, liion, h2tank, elyz, fc, tes, heater, grid, ω_simu)

            # Compute operation decision variables
            compute_operation_decisions(h, y, s, ld, pv, liion, h2tank, elyz, fc, tes, heater, controller, ω_optim, parameters)

            # Compute operation dynamics
            compute_operation_dynamics(h, y, s, liion, h2tank, elyz, fc, tes, heater, controller, parameters)

            # Power balance constraint checked for each node
            power_balance(h, y, s, ld, pv, liion, h2tank, elyz, fc, tes, heater, grid)

        end

        # Update investment informations
        update_investment_informations(y, s, pv, liion, h2tank, elyz, fc, tes, heater, ω_simu)

        # Compute investment decision variables
        compute_investment_decisions(y, s, ld, pv, liion, h2tank, elyz, fc, tes, heater, designer, ω_optim, parameters)

        # Compute investment dynamics
        compute_investment_dynamics(y, s, pv, liion, h2tank, elyz, fc, tes, designer)

    end
end
