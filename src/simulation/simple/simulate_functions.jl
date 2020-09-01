#=
    This file includes all the funtions needed to simulate the DES
 =#

# Simulate function
function simulate(ld::Load, pv::Source, liion::Liion, controller::AbstractController,
     designer::AbstractDesigner, grid::Grid, ω_optim::Scenarios, ω_simu::Scenarios, parameters::NamedTuple)

    # Parameters
    ns = size(ld.power_E,3)

    # We simulate over the horizon for all the scenarios
    @showprogress for s in 1:ns
        simulate_scenario(s, ld, pv, liion, controller,
             designer, grid, ω_optim, ω_simu, parameters)
    end

    # Compute economic indicators
    costs = compute_economics(ld, pv, liion, designer, grid, parameters)

    return costs
end

# Parallel simulation
function simulate_parallel(ld::Load, pv::Source, liion::Liion, controller::AbstractController,
     designer::AbstractDesigner, grid::Grid, ω_optim::Scenarios, ω_simu::Scenarios, parameters::NamedTuple)

    # Parameters
    ns = size(ld.power_E,3)

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

    # Compute economic indicators
    costs = compute_economics(ld, pv, liion, designer, grid, parameters)

    return costs
end

# Simulate function
function simulate_parallel_distributed(ld::Load, pv::Source, liion::Liion, controller::AbstractController,
     designer::AbstractDesigner, grid::Grid, ω_optim::Scenarios, ω_simu::Scenarios, parameters::NamedTuple)

    # Parameters
    ns = size(ld.power_E,3)

    # We simulate over the horizon for all the scenarios in parallel using the distributed macro
    @sync @distributed for s in 1:ns
        simulate_scenario(s, ld, pv, liion, controller,
             designer, grid, ω_optim, ω_simu, parameters)
    end

    # Compute economic indicators
    costs = compute_economics(ld, pv, liion, designer, grid, parameters)

    return costs
end

# Simulate function thread
function simulate_parallel_thread(ld::Load, pv::Source, liion::Liion, controller::AbstractController,
     designer::AbstractDesigner, grid::Grid, ω_optim::Scenarios, ω_simu::Scenarios, parameters::NamedTuple)

    # Parameters
    ns = size(ld.power_E,3)

    # We simulate over the horizon for all the scenarios in parallel using the distributed macro
    Threads.@threads for s in 1:ns
        simulate_scenario(s, ld, pv, liion, controller,
             designer, grid, ω_optim, ω_simu, parameters)
    end

    # Compute economic indicators
    costs = compute_economics(ld, pv, liion, designer, grid, parameters)

    return costs
end

# Simulate one scenario
function simulate_scenario(s::Int64, ld::Load, pv::Source, liion::Liion, controller::AbstractController,
     designer::AbstractDesigner, grid::Grid, ω_optim::Scenarios, ω_simu::Scenarios, parameters::NamedTuple)

    # Parameters
    nh = size(ld.power_E,1)
    ny = size(ld.power_E,2)

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
                power_balance_checking(h, y, s, ld, pv, liion, grid)

            end

            # Update investment informations
            update_investment_informations(y, s, pv, liion, ω_simu)

            # Compute investment decision variables
            compute_investment_decisions(y, s, ld, pv, liion, grid, controller, designer, ω_optim, parameters)

            # Compute investment dynamics
            compute_investment_dynamics(y, s, pv, liion, designer)

    end

end
