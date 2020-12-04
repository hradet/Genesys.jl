#=
    Designer based on a metaheuristic
=#

mutable struct MetaheuristicOptions
    method::Metaheuristics.AbstractMetaheuristic
    iterations::Int64
    controller::AbstractController
    isnpv::Bool
    risk_measure::String
    reducer::AbstractScenariosReducer
    share_constraint::String
    lpsp_constraint::String
    tol_lpsp::Float64
    reopt::Bool

    MetaheuristicOptions(; method = Metaheuristics.Clearing(),
                           iterations = 20,
                           controller = RBC(),
                           isnpv = false,
                           risk_measure = "esperance",
                           reducer = KmeansReducer(),
                           share_constraint = "hard",
                           lpsp_constraint = "soft",
                           tol_lpsp = 0.05,
                           reopt = false) =
                           new(method, iterations, controller, isnpv, risk_measure, reducer, share_constraint, lpsp_constraint, tol_lpsp, reopt)

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
    nh = size(ω.ld_E.power,1)
    ny = size(ω.ld_E.power,2)
    ns = size(ω.ld_E.power,3)

    # Initialize DES
    des_m = copy(des, nh, ny, ns)

    # Initialize controller
    controller_m = initialize_controller!(des_m, designer.options.controller, ω)

    # Initialize with the dummy designer
    designer_m = initialize_designer!(des_m, DummyDesigner(), ω)

    # Initialize with the decisions variables
    designer_m.u.pv[1,:] .= decisions[1]
    designer_m.u.liion[1,:] .= decisions[2]
    designer_m.u.h2tank[1,:] .= decisions[3]
    designer_m.u.elyz[1,:] .= decisions[4]
    designer_m.u.fc[1,:] .= decisions[5]
    designer_m.u.tes[1,:] .= decisions[6]

    # Simulate
    simulate!(des_m, controller_m, designer_m, ω, options = Genesys.Options(mode = "multithreads"))

    # Objective - algorithm find the maximum
    if designer.options.isnpv
        obj = sum(probabilities[s] * Costs(des_m,designer_m).npv[s] for s in 1:ns)
    else
        obj = - compute_annualised_capex(1, 1, des_m, designer_m) - sum(probabilities[s] * compute_grid_cost(2, s, des_m) for s in 1:ns)
    end

    # LPSP constraint for the heat
    if isa(des_m.ld_H, Load)
        if designer.options.lpsp_constraint == "soft"
            obj -= 1e32 * max(0., sum(probabilities[s] * LPSP(y, s, des_m).lpsp_H - designer.options.tol_lpsp for s in 1:ns, y in 2:ny))
        elseif designer.options.lpsp_constraint == "hard"
            obj -= 1e32 * max(0., maximum(LPSP(y, s, des_m).lpsp_H - designer.options.tol_lpsp for s in 1:ns, y in 2:ny))
        end
    end

    # SoC constraint for the seasonal storage
    isa(des_m.h2tank, H2Tank) ? obj -= 1e32 *  max(0., maximum(des_m.h2tank.soc[1,y,s] - des_m.h2tank.soc[end,y,s] for s in 1:ns, y in 2:ny)) : nothing

    # Share constraint
    if designer.options.share_constraint == "hard"
        obj -= 1e32 * max(0., des.parameters.τ_share - minimum(compute_share(y, s, des_m) for s in 1:ns, y in 2:ny))
    elseif designer.options.share_constraint == "soft"
        obj -= 1e32 * max(0., des.parameters.τ_share - sum(probabilities[s] * compute_share(y, s, des_m) for s in 1:ns, y in 2:ny))
    end

    return obj
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
        ω_reduced, probabilities = reduce(designer.options.reducer, ω)
        # Repeat to simulate 2 years
        ω_reduced = repeat(ω_reduced, 1, 2, 1)
    end

    # Bounds
    lb, ub = set_bounds(des)

    # Optimize
    designer.results = Metaheuristics.optimize(lb, ub,
                                               designer.options.method,
                                               options = Metaheuristics.Options(iterations=designer.options.iterations, multithreads=true)
    ) do decisions
        fobj(decisions, des, designer, ω_reduced, probabilities)
      end

    # Assign values
    designer.u.pv[1,:] .= designer.results.minimizer[1]
    designer.u.liion[1,:] .= designer.results.minimizer[2]
    designer.u.h2tank[1,:] .= designer.results.minimizer[3]
    designer.u.elyz[1,:] .= designer.results.minimizer[4]
    designer.u.fc[1,:] .= designer.results.minimizer[5]
    designer.u.tes[1,:] .= designer.results.minimizer[6]

    # Save history for online optimization
    designer.history = ω_reduced

    return designer
end

### Online
function compute_investment_decisions!(y::Int64, s::Int64, des::DistributedEnergySystem, designer::Metaheuristic)
    ϵ = 0.1

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

    lb, ub = zeros(6), zeros(6)
    isa(des.pv, Source) ? ub[1] = 1000. : nothing
    isa(des.liion, Liion) ? ub[2] = 1000. : nothing
    isa(des.h2tank, H2Tank) ? ub[3] = 50000. : nothing
    isa(des.elyz, Electrolyzer) ? ub[4] = 50. : nothing
    isa(des.fc, FuelCell) ? ub[5] = 50. : nothing
    isa(des.tes, ThermalSto) ? ub[6] = 1000. : nothing

    return lb, ub
end
