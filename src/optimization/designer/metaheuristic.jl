#=
    Designer based on a metaheuristic
=#

mutable struct MetaheuristicOptions
    method::Metaheuristics.AbstractMetaheuristic
    iterations::Int64
    multithreads::Bool
    controller::AbstractController
    isnpv::Bool
    reducer::AbstractScenariosReducer
    objective_risk::AbstractRiskMeasure
    share_risk::AbstractRiskMeasure
    lpsp_risk::AbstractRiskMeasure
    lpsp_tol::Float64
    reopt::Bool
    read_reduction::Union{String, Nothing}
    write_reduction::Union{String, Nothing}

    MetaheuristicOptions(; method = Metaheuristics.Clearing(),
                           iterations = 100,
                           multithreads = false,
                           controller = RBC(),
                           isnpv = false,
                           reducer = FeatureBasedReducer(),
                           objective_risk = Expectation(),
                           share_risk = Expectation(),
                           lpsp_risk = WorstCase(),
                           lpsp_tol = 1e-3,
                           reopt = false,
                           read_reduction = nothing,
                           write_reduction = nothing) =
                           new(method, iterations, multithreads, controller, isnpv, reducer, objective_risk, share_risk, lpsp_risk, lpsp_tol, reopt, read_reduction, write_reduction)

end

mutable struct Metaheuristic <: AbstractDesigner
    options::MetaheuristicOptions
    u::NamedTuple
    results::Metaheuristics.MetaheuristicResults
    history::AbstractScenarios

    Metaheuristic(; options = MetaheuristicOptions()) = new(options)
end

# Objective functions
function fobj(decisions::Array{Float64,1}, des::DistributedEnergySystem, designer::Metaheuristic, ω::Scenarios, probabilities::Array{Float64})
    # Paramters
    nh, ny, ns = size(ω.ld_E.power)
    λ1 = λ2 = λ3 = 1e6

    # Initialize DES
    des_m = copy(des, nh, ny, ns)

    # Initialize controller
    controller_m = initialize_controller!(des_m, designer.options.controller, ω)

    # Initialize with the manual designer
    designer_m = initialize_designer!(des_m, Manual(pv = decisions[1], liion = decisions[2], tes = decisions[3], h2tank = decisions[4], elyz = decisions[5], fc = decisions[6]), ω)

    # Simulate
    simulate!(des_m, controller_m, designer_m, ω)

    # Metrics
    metrics = Metrics(des_m, designer_m)

    # Share constraint
    share = max(0., des.parameters.renewable_share - conditional_value_at_risk([reshape(metrics.renewable_share[2:ny, 1:ns], :, 1)...],  probabilities,  designer.options.share_risk))

    # LPSP constraint for the heat
    isa(des_m.ld_H, Load) ? lpsp = max(0., conditional_value_at_risk([reshape(metrics.lpsp.heat[2:ny, 1:ns], :, 1)...],  probabilities,  designer.options.lpsp_risk) - designer.options.lpsp_tol) : lpsp = 0.

    # SoC constraint for the seasonal storage
    isa(des_m.h2tank, H2Tank) ? soc_h2tank = sum(max(0., des_m.h2tank.soc[1,y,s] - des_m.h2tank.soc[end,y,s]) for y in 2:ny, s in 1:ns) : soc_h2tank = 0.

    # Objective - Algortihm find the maximum
    if designer.options.isnpv
        # NPV
        npv = conditional_value_at_risk([metrics.npv.total...], probabilities, designer.options.objective_risk)
        return npv - λ1 * share - λ2 * lpsp - λ3 * soc_h2tank
    else
        # Equivalent annual cost
        eac = conditional_value_at_risk([metrics.eac.total...], probabilities, designer.options.objective_risk)
        return - eac - λ1 * share - λ2 * lpsp - λ3 * soc_h2tank
    end
end

### Offline
function initialize_designer!(des::DistributedEnergySystem, designer::Metaheuristic, ω::Scenarios{Array{DateTime,3}, Array{Float64,3}, Array{Float64,2}})
    # Preallocate and assigned values
    preallocate!(designer, des.parameters.ny, des.parameters.ns)

    # Scenario reduction from the optimization scenario pool
    println("Starting scenario reduction...")
    if designer.options.isnpv
        ω_reduced, probabilities = reduce(designer.options.reducer, ω)
    else
        if isa(designer.options.read_reduction, Nothing)
            println("Starting scenario reduction...")
            ω_reduced, probabilities = reduce(designer.options.reducer, ω)
            # Saving
            if !isa(designer.options.write_reduction, Nothing)
                save(designer.options.write_reduction, "scenarios", ω_reduced, "probabilities", probabilities)
            end
        else
            println("Reading scenario reduction from file...")
            ω_reduced = load(designer.options.read_reduction, "scenarios")
            probabilities = load(designer.options.read_reduction, "probabilities")
        end
        # Repeat to simulate 2 years
        ω_reduced = repeat(ω_reduced, 1, 2, 1)
    end

    # Bounds
    lb, ub = set_bounds(des)

    # Optimize
    designer.results = Metaheuristics.optimize(lb, ub,
                                               designer.options.method,
                                               options = Metaheuristics.Options(iterations = designer.options.iterations, multithreads = designer.options.multithreads)
    ) do decisions
        fobj(decisions, des, designer, ω_reduced, probabilities)
      end

    # Assign values
    designer.u.pv[1,:] .= designer.results.minimizer[1]
    designer.u.liion[1,:] .= designer.results.minimizer[2]
    designer.u.tes[1,:] .= designer.results.minimizer[3]
    designer.u.h2tank[1,:] .= designer.results.minimizer[4]
    designer.u.elyz[1,:] .= designer.results.minimizer[5]
    designer.u.fc[1,:] .= designer.results.minimizer[6]

    # Save history for online optimization
    designer.history = ω_reduced

    return designer
end

### Online
function compute_investment_decisions!(y::Int64, s::Int64, des::DistributedEnergySystem, designer::Metaheuristic)
    ϵ = 0.2

    if designer.options.reopt && y != 1
        # Do we need to reoptimize ?
        (isa(des.liion, Liion) && des.liion.soh[end,y,s] < ϵ) || (isa(des.elyz, Electrolyzer) && des.elyz.soh[end,y,s] < ϵ) || (isa(des.fc, FuelCell) && des.fc.soh[end,y,s] < ϵ) ? nothing : return
        println("Re-optimization not yet implemented...")
    else
        isa(des.liion, Liion) && des.liion.soh[end,y,s] < ϵ ? designer.u.liion[y,s] = designer.u.liion[1,s] : nothing
        isa(des.elyz, Electrolyzer) && des.elyz.soh[end,y,s] < ϵ ? designer.u.elyz[y,s] = designer.u.elyz[1,s] : nothing
        isa(des.fc, FuelCell) && des.fc.soh[end,y,s] < ϵ ? designer.u.fc[y,s] = designer.u.fc[1,s] : nothing
    end
end

### Utils
function set_bounds(des::DistributedEnergySystem)
    lb, ub = 1e-6 * ones(6), zeros(6)
    isa(des.pv, Source) ? ub[1] = 1000. : nothing
    isa(des.liion, Liion) ? ub[2] = 1000. : nothing
    isa(des.tes, ThermalSto) ? ub[3] = 1000. : nothing
    isa(des.h2tank, H2Tank) ? ub[4] = 50000. : nothing
    isa(des.elyz, Electrolyzer) ? ub[5] = 50. : nothing
    isa(des.fc, FuelCell) ? ub[6] = 50. : nothing
    return lb, ub
end
