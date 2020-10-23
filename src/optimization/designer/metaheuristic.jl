#=
    Designer based on a metaheuristic
=#

mutable struct MetaheuristicOptions
    method::Metaheuristics.AbstractMetaheuristic
    controller::AbstractController
    iterations::Int64
    scenario_reduction::String
    share_constraint::Bool
    reopt::Bool
    s::Int64

    MetaheuristicOptions(; method = Metaheuristics.Clearing(),
                           controller = RBC(),
                           iterations = 10,
                           scenario_reduction = "manual",
                           share_constraint = true,
                           reopt = false,
                           s = 1) =
                           new(method, controller, iterations, scenario_reduction, share_constraint, reopt, s)

end

mutable struct Metaheuristic <: AbstractMultiStageDesigner
    options::MetaheuristicOptions
    u::NamedTuple
    results::Metaheuristics.MetaheuristicResults
    history::AbstractScenarios

    Metaheuristic(; options = MetaheuristicOptions()) = new(options)
end

# Objective functions
function fobj(decisions, des, designer, ω)
    # Initialize DES
    des_m = deepcopy(des)

    # Initialize controller
    controller_m = initialize_controller!(des_m, designer.options.controller, ω)

    # Initialize with the dummy designer
    designer_m = initialize_designer!(des_m, DummyDesigner(), ω)

    # Initialize with the decisions variables
    #TODO bug si il en manque un dans le lot... decision 4 devient decision 3...
    isa(des.pv, Source) ? designer_m.u.pv[1,:] .= decisions[1] : nothing
    isa(des.liion, Liion) ? designer_m.u.liion[1,:] .= decisions[2]  : nothing
    isa(des.h2tank, H2Tank) ? designer_m.u.h2tank[1,:] .= decisions[3]  : nothing
    isa(des.elyz, Electrolyzer) ? designer_m.u.elyz[1,:] .= decisions[4]  : nothing
    isa(des.fc, FuelCell) ? designer_m.u.fc[1,:] .= decisions[5]  : nothing
    isa(des.tes, ThermalSto) ? designer_m.u.tes[1,:] .= decisions[6]  : nothing

    # Simulate
    simulate!(1, des_m, controller_m, designer_m, ω)

    # Compute metrics
    metrics_m = compute_metrics(1, des_m, designer_m)

    # Objective - algorithm find the maximum
   designer.options.share_constraint ? obj = metrics_m.costs.npv[1] - 1e32 * max(0., des.parameters.τ_share - minimum(metrics_m.τ_share[2:end,:])) : obj = metrics_m.costs.npv[1]

    return obj
end

### Offline
function initialize_designer!(des::DistributedEnergySystem, designer::Metaheuristic, ω_optim::AbstractScenarios)

    # Save history for online optimization
    designer.history = ω_optim

    # Preallocate and assigned values
    preallocate!(designer, des.parameters.ny, des.parameters.ns)

    return designer
end

### Online
function compute_investment_decisions!(y::Int64, s::Int64, des::DistributedEnergySystem, designer::Metaheuristic)
    ϵ = 0.1

    if s == 1 && y == 1
        # Scenario reduction from the optimization scenario pool
        ω_meta = Genesys.scenarios_reduction(designer, designer.history)

        # Bounds
        lb, ub = Genesys.set_bounds(des)

        # Optimize
        designer.results = Metaheuristics.optimize(lb, ub,
                                                   designer.options.method,
                                                   options = Metaheuristics.Options(iterations=designer.options.iterations, multithreads=true)
        ) do decisions
            fobj(decisions, des, designer, ω_meta)
          end

        # Assign values
        isa(des.pv, Source) ? designer.u.pv[1,:] .= designer.results.minimizer[1] : nothing
        isa(des.liion, Liion) ? designer.u.liion[1,:] .= designer.results.minimizer[2] : nothing
        isa(des.h2tank, H2Tank) ? designer.u.h2tank[1,:] .= designer.results.minimizer[3] : nothing
        isa(des.elyz, Electrolyzer) ? designer.u.elyz[1,:] .= designer.results.minimizer[4] : nothing
        isa(des.fc, FuelCell) ? designer.u.fc[1,:] .= designer.results.minimizer[5] : nothing
        isa(des.tes, ThermalSto) ? designer.u.tes[1,:] .= designer.results.minimizer[5] : nothing

    elseif designer.options.reopt
        # Do we need to reoptimize ?
        (isa(des.liion, Liion) && des.liion.soh[end,y,s] < ϵ) || (isa(des.elyz, Electrolyzer) && des.elyz.soh[end,y,s] < ϵ) || (isa(des.fc, FuelCell) && des.fc.soh[end,y,s] < ϵ) ? nothing : return
        println("Re-optimization not yet implemented...")
    else
        isa(des.liion, Liion) && des.liion.soh[end,y,s] < ϵ ? designer.u.liion[y,s] = designer.u.liion[1,s] : nothing
        isa(des.elyz, Electrolyzer) && des.elyz.soh[end,y,s] < ϵ ? designer.u.elyz[y,s] = designer.u.elyz[1,s] : nothing
        isa(des.fc, FuelCell) && des.fc.soh[end,y,s] < ϵ ? designer.u.fc[y,s] = designer.u.fc[1,s] : nothing
    end
end
