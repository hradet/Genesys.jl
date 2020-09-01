#=
    This file includes all the funtions needed to simulate the DES
 =#

# Simulate function
function simulate(ld::Load, pv::Source, liion::Liion, h2tank::H2Tank,
    elyz::Electrolyzer, fc::FuelCell, tes::ThermalSto, heater::Heater,
    controller::AbstractController, designer::AbstractDesigner, grid::Grid,
    ω_optim::Scenarios, ω_simu::Scenarios, parameters::NamedTuple)

    # Parameters
    nh = size(ld.power_E,1)
    ny = size(ld.power_E,2)
    ns = size(ld.power_E,3)

    # We simulate over the horizon for all the scenarios
    for s in 1:ns
        for y = 1 : ny
                for h = 1 : nh

                    # Update operation informations
                    update_operation_informations(h, y, s, ld, pv, liion, h2tank, elyz, fc, tes, heater, grid, ω_simu)

                    # Compute operation decision variables
                    compute_operation_decisions(h, y, s, ld, pv, liion, h2tank, elyz, fc, tes, heater, controller, ω_optim, parameters)

                    # Compute operation dynamics
                    compute_operation_dynamics(h, y, s, liion, h2tank, elyz, fc, tes, heater, controller, parameters)

                    # Power balance constraint checked for each node
                    power_balance_checking(h, y, s, ld, pv, liion, h2tank, elyz, fc, tes, heater, grid)

                end

                # Update investment informations
                update_investment_informations(y, s, pv, liion, h2tank, elyz, fc, tes, heater, ω_simu)

                # Compute investment decision variables
                compute_investment_decisions(y, s, ld, pv, liion, h2tank, elyz, fc, tes, heater, designer, ω_optim, parameters)

                # Compute investment dynamics
                compute_investment_dynamics(y, s, pv, liion, h2tank, elyz, fc, tes, designer)

        end
    end

    # Compute economic indicators
    costs = compute_economics(ld, pv, liion, h2tank, elyz, fc, tes, heater, designer, grid, parameters)

    return costs
end

# Parallel simulation
# TODO
