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
    obj::String
    s::Int64

    MetaheuristicOptions(; method = Metaheuristics.Clearing(),
                           controller = RBC(),
                           iterations = 20,
                           scenario_reduction = "manual",
                           share_constraint = true,
                           reopt = false,
                           obj = "npv",
                           s = 1) =
                           new(method, controller, iterations, scenario_reduction, share_constraint, reopt, obj, s)

end

mutable struct Metaheuristic <: AbstractMultiStageDesigner
    options::MetaheuristicOptions
    u::NamedTuple
    results::Metaheuristics.MetaheuristicResults
    history::AbstractScenarios

    Metaheuristic(; options = MetaheuristicOptions()) = new(options)
end

# Objective functions
function fobj_npv(decisions, des, designer, ω_m)
    # Initialize DES
    des_m = deepcopy(des)

    # Initialize controller
    controller_m = initialize_controller!(des_m, designer.options.controller, ω_m)

    # Initialize with the dummy designer
    designer_m = initialize_designer!(des_m, DummyDesigner(), ω_m)

    # Initialize with the decisions variables
    designer_m.u.pv[1,:] .= decisions[1]
    designer_m.u.liion[1,:] .= decisions[2]
    designer_m.u.h2tank[1,:] .= decisions[3]
    designer_m.u.elyz[1,:] .= decisions[4]
    designer_m.u.fc[1,:] .= decisions[5]
    designer_m.u.tes[1,:] .= decisions[6]

    # Simulate
    simulate!(1, des_m, controller_m, designer_m, ω_m, Options())

    # Metrics
    metrics_m = Metrics(1, des_m, designer_m)

    # Objective - algorithm find the maximum
   obj = metrics_m.costs.npv[1]

   # Add the LPSP constraint for the heat
   isa(des.ld_H, Load) ? obj -= 1e32 * max(0., maximum(metrics_m.lpsp.lpsp_H[2:end,:]) - 0.05) : nothing

   # Add the share constraint
   designer.options.share_constraint ? obj -= 1e32 * max(0., des.parameters.τ_share - minimum(metrics_m.τ_share[2:end,:])) : nothing

    return obj
end
function fobj_eac(decisions, des, designer, ω_m)
    # Initialize DES
    des_m = deepcopy(des)

    # Initialize controller
    controller_m = initialize_controller!(des_m, designer.options.controller, ω_m)

    # Initialize with the dummy designer
    designer_m = initialize_designer!(des_m, DummyDesigner(), ω_m)

    # Initialize with the decisions variables
    designer_m.u.pv[1,:] .= decisions[1]
    designer_m.u.liion[1,:] .= decisions[2]
    designer_m.u.h2tank[1,:] .= decisions[3]
    designer_m.u.elyz[1,:] .= decisions[4]
    designer_m.u.fc[1,:] .= decisions[5]
    designer_m.u.tes[1,:] .= decisions[6]

    # Simulate
    for y in 1:2
        simulate!(y, 1, des_m, controller_m, designer_m, ω_m, Options())
    end

    # Objective - algorithm find the maximum
   obj = - compute_annualised_capex(1, 1, des_m, designer_m) - compute_grid_cost(2, 1, des_m)

   # Add the LPSP constraint for the heat
   isa(des_m.ld_H, Load) ? obj -= 1e32 * max(0., Genesys.LPSP(2, 1, des_m).lpsp_H - 0.05) : nothing

   # Add the share constraint
   designer.options.share_constraint ? obj -= 1e32 * max(0., des.parameters.τ_share - compute_share(2, 1, des_m)) : nothing

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
        ω_meta = scenarios_reduction(designer, designer.history)

        # Bounds
        lb, ub = set_bounds(des)

        # Optimize
        designer.results = Metaheuristics.optimize(lb, ub,
                                                   designer.options.method,
                                                   options = Metaheuristics.Options(iterations=designer.options.iterations, multithreads=true)
        ) do decisions
            if designer.options.obj == "npv"
                fobj_npv(decisions, des, designer, ω_meta)
            elseif designer.options.obj == "eac"
                fobj_eac(decisions, des, designer, ω_meta)
            else
                println("Objective function unknown...")
            end
          end

        # Assign values
        designer.u.pv[1,:] .= designer.results.minimizer[1]
        designer.u.liion[1,:] .= designer.results.minimizer[2]
        designer.u.h2tank[1,:] .= designer.results.minimizer[3]
        designer.u.elyz[1,:] .= designer.results.minimizer[4]
        designer.u.fc[1,:] .= designer.results.minimizer[5]
        designer.u.tes[1,:] .= designer.results.minimizer[6]

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
