#=
    This file includes all the funtions needed to simulate the DES
 =#

 mutable struct Options
     mode::String
     Options(; mode = "serial") = new(mode)
 end

# Simulate function
function simulate!(des::DistributedEnergySystem,
                   controller::AbstractController,
                   designer::AbstractDesigner,
                   ω_simu::AbstractScenarios;
                   options = Options())

    # Parameters
    ns = des.parameters.ns

    if options.mode == "serial"
        # We simulate over the horizon for all the scenarios
        @showprogress for s in 1:ns
            simulate_scenario!(s, des, controller, designer, ω_simu)
        end
    elseif options.mode == "multicores"
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
                        remotecall_fetch(simulate_scenario!, p, s, des, controller, designer, ω_simu)
                    end
                end
            end
        end
    elseif options.mode == "distributed"
        # We simulate over the horizon for all the scenarios in parallel using the distributed macro
        @sync @distributed for s in 1:ns
            simulate_scenario!(s, des, controller, designer, ω_simu)
        end
    elseif options.mode == "multithreads"
        # We simulate over the horizon for all the scenarios in parallel using the distributed macro
        Threads.@threads for s in 1:ns
            simulate_scenario!(s, des, controller, designer, ω_simu)
        end
    else
        println("Unknown mode... Please chose between 'serial', 'multicores', 'distributed' or 'multithreads'.")
    end
end

function simulate_scenario!(s::Int64,
                            des::DistributedEnergySystem,
                            controller::AbstractController,
                            designer::AbstractDesigner,
                            ω_simu::AbstractScenarios)

    # Parameters
    nh = des.parameters.nh
    ny = des.parameters.ny

    # We simulate over the horizon for a signle scenario
    for y = 1 : ny
        for h = 1 : nh

            # Update operation informations
            update_operation_informations!(h, y, s, des, ω_simu)

            # Compute operation decision variables
            compute_operation_decisions!(h, y, s, des, controller)

            # Compute operation dynamics
            compute_operation_dynamics!(h, y, s, des, controller)

            # Power balance constraint checked for each node
            compute_power_balances!(h, y, s, des)

        end

        # Update investment informations
        update_investment_informations!(y, s, des, ω_simu)

        # Compute investment decision variables
        compute_investment_decisions!(y, s, des, designer)

        # Compute investment dynamics
        compute_investment_dynamics!(y, s, des, designer)

    end
end
