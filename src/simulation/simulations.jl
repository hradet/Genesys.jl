#=
    This file includes all the funtions needed to simulate the DES
 =#

 mutable struct Options
     mode::String
     Options(; mode = "serial") = new(mode)
 end

# Simulate function
function simulate!(des::DES, ω_simu::AbstractScenarios; options = Options())

    # Parameters
    ns = size(ω_simu.values.ld_E,3)

    if options.mode == "serial"
        # We simulate over the horizon for all the scenarios
        @showprogress for s in 1:ns
            simulate_scenario!(s, des, ω_simu)
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
                        remotecall_fetch(simulate_scenario!, p, s, des, ω_simu)
                    end
                end
            end
        end
    elseif options.mode == "distributed"
        # We simulate over the horizon for all the scenarios in parallel using the distributed macro
        @sync @distributed for s in 1:ns
            simulate_scenario!(s, des, ω_simu)
        end
    elseif options.mode == "multithreads"
        # We simulate over the horizon for all the scenarios in parallel using the distributed macro
        Threads.@threads for s in 1:ns
            simulate_scenario!(s, des, ω_simu)
        end
    else
        println("Unknown mode... Please chose between 'serial', 'multicores', 'distributed' or 'multithreads'.")
    end
end

function simulate_scenario!(s::Int64, des::DES, ω_simu::AbstractScenarios)

    # Parameters
    nh = size(ω_simu.values.ld_E,1)
    ny = size(ω_simu.values.ld_E,2)

    # We simulate over the horizon for a signle scenario
    for y = 1 : ny
        for h = 1 : nh

            # Update operation informations
            update_operation_informations!(h, y, s, des, ω_simu)

            # Compute operation decision variables
            compute_operation_decisions!(h, y, s, des)

            # Compute operation dynamics
            compute_operation_dynamics!(h, y, s, des)

            # Power balance constraint checked for each node
            compute_power_balances!(h, y, s, des)

        end

        # Update investment informations
        update_investment_informations!(y, s, des, ω_simu)

        # Compute investment decision variables
        compute_investment_decisions!(y, s, des)

        # Compute investment dynamics
        compute_investment_dynamics!(y, s, des)

    end
end
