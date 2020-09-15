#=
    Designer based on a metaheuristic
=#
#TODO ajouter controller dans les options de metheuristic
mutable struct MetaHeuristicDesigner <: AbstractMultiStageDesigner
    u::NamedTuple
    parameters::Dict{String, Any}
    results
    MetaHeuristicDesigner() = new()
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

     # Define the objective function from global variables
     function objective_function(u::Array{Float64,1})

          # Initialize controller
          dummy_designer = DummyDesigner()

          # Initialize designer
          initialize_designer(ld, pv, liion, dummy_designer, grid, ω_meta, parameters)
          dummy_designer.u.u_pv[1], dummy_designer.u.u_liion[1],  = u[1], u[2]

          # Simulate
          simulate(ld, pv, liion, controller, dummy_designer, grid, ω_meta, ω_meta, parameters)

          # Compute cost indicators
          costs = compute_economics(ld, pv, liion, dummy_designer, grid, parameters)

          # Compute tech indicators
          tech = compute_tech_indicators(ld, grid)

          # Objective - algorithm find the maximum
          obj = costs.npv[1] - 1e10 * max(0., grid.τ_energy - minimum(tech.τ_self[2:end,:]))

         return obj
     end

     # Initialize decision variables and bounds
     lb = [designer.parameters["lb"].pv, designer.parameters["lb"].liion]
     ub = [designer.parameters["ub"].pv, designer.parameters["ub"].liion]

     # Compute investment decisions
     designer.results = Metaheuristics.optimize(objective_function, lb, ub,
                                                Metaheuristics.Clearing(),
                                                options = Metaheuristics.Options(iterations=10, multithreads=true, itmax_unchanged=5))

     # Formatting variables to simulation
     designer.u = (
     u_pv = repeat(vcat(designer.results.minimizer[1], zeros(ny-1,1)), 1, ns),
     u_liion = repeat(vcat(designer.results.minimizer[2], zeros(ny-1,1)), 1, ns),
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
