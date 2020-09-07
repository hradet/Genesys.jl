#=
    Designer based on a metaheuristic
=#

mutable struct MetaHeuristicDesigner <: AbstractMultiStageDesigner
    u::NamedTuple
    parameters::Dict{String, Any}
    results
    MetaHeuristicDesigner() = new()
end

### Models ###
# Simple
function simulation_metaheuristic(u::Array{Float64,1})

     # GUI loading
     outputGUI = loadGUI("RuleBasedController", "DummyDesigner", ns=1)

     # Initialization without any controller and designer
     ld, pv, liion, _, _, _, _, _, controller, designer, grid, ω_optim, ω_simu = initialization(outputGUI)

     # Initialize controller
     initialize_controller(ld, pv, liion, controller, grid, ω_optim, outputGUI["parameters"])

     # Initialize designer
     initialize_designer(ld, pv, liion, designer, grid, ω_optim, outputGUI["parameters"])
     designer.u.u_pv[1], designer.u.u_liion[1],  = u[1], u[2]

     # Simulate
     simulate(ld, pv, liion, controller, designer, grid, ω_optim, ω_simu, outputGUI["parameters"])

    costs = compute_economics(ld, pv, liion, designer, grid, outputGUI["parameters"])

    return -costs.npv[1] # As the algorithm find the minimum...
end

#### Offline functions ####
# Simple
function initialize_designer(ld::Load, pv::Source, liion::Liion,
     controller::AbstractController, designer::MetaHeuristicDesigner, grid::Grid, ω_optim::Scenarios,
     parameters::NamedTuple)

     # Parameters
     ny = size(ld.power_E,2) # number of simulation years
     ns = size(ld.power_E,3) # number of scenarios

     # Scenario reduction from the optimization scenario pool
     ω_meta = scenarios_reduction(designer, ω_optim)

     # Initialize decision variables
     u_init = [designer.parameters["u_init"].u_pv, designer.parameters["u_init"].u_liion]

     # Compute investment decisions
     designer.results = Evolutionary.optimize(simulation_metaheuristic, [0., 0.], [500., 500.], CMAES(mu=20, lambda=100))

     # Formatting variables to simulation
     designer.u = (
     u_pv = repeat(vcat(Evolutionary.minimizer(results)[1], zeros(ny-1,1)), 1, ns),
     u_liion = repeat(vcat(Evolutionary.minimizer(results)[2], zeros(ny-1,1)), 1, ns),
     )

end

#### Online functions ####
# Simple
function compute_investment_decisions(y::Int64, s::Int64, ld::Load, pv::Source,
    liion::Liion, grid::Grid, designer::MetaHeuristicDesigner, ω_optim::Scenarios, parameters::NamedTuple)
    ϵ = 0.1
    if liion.soh[end,y,s] < ϵ
        designer.u.u_liion[y,s] = liion.Erated[y,s]
    end
end
# Multi-energy
function compute_investment_decisions(y::Int64, s::Int64, ld::Load, pv::Source,
    liion::Liion, h2tank::H2Tank, elyz::Electrolyzer, fc::FuelCell, tes::ThermalSto,
    heater::Heater, designer::MetaHeuristicDesigner, ω_optim::Scenarios, parameters::NamedTuple)
    ϵ = 0.1

    # Liion
    if liion.soh[end,y,s] < ϵ
        designer.u.u_liion[y,s] = liion.Erated[y,s]
    end

    # Electrolyzer
    if elyz.soh[end,y,s] < ϵ
        designer.u.u_elyz[y,s] = elyz.powerMax[y,s]
    end

    # FuelCell
    if fc.soh[end,y,s] < ϵ
        designer.u.u_fc[y,s] = fc.powerMax[y,s]
    end
end
