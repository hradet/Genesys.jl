#=
    Designer based on a metaheuristic
=#

mutable struct MetaheuristicOptions
    method::Metaheuristics.AbstractMetaheuristic
    iterations::Int64
    controller::AbstractController
    mode::String
    risk_measure::String
    scenario_reduction::Bool
    share_constraint::Bool
    reopt::Bool

    MetaheuristicOptions(; method = Metaheuristics.Clearing(),
                           iterations = 20,
                           controller = RBC(),
                           mode = "deterministic", # "deterministic", "twostage" or "npv"
                           risk_measure = "esperance",
                           scenario_reduction = true,
                           share_constraint = true,
                           reopt = false) =
                           new(method, iterations, controller, mode, risk_measure, scenario_reduction, share_constraint, reopt)

end

mutable struct Metaheuristic <: AbstractDesigner
    options::MetaheuristicOptions
    u::NamedTuple
    results::Metaheuristics.MetaheuristicResults
    history::AbstractScenarios

    Metaheuristic(; options = MetaheuristicOptions()) = new(options)
end

# Objective functions
function fobj(decisions::Array{Float64,1}, des::DistributedEnergySystem, designer::Metaheuristic, ω_m::AbstractScenarios)
    # Paramters
    nh = size(ω_m.ld_E.power,1)
    ny = size(ω_m.ld_E.power,2)
    ns = size(ω_m.ld_E.power,3)

    # Initialize DES
    des_m = copy(des, nh, ny, ns)

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
    for s in 1:ns
        simulate!(s, des_m, controller_m, designer_m, ω_m, Genesys.Options())
    end

    # Objective - algorithm find the maximum
    if designer.options.mode == "npv"
        obj = mean(Costs(des_m,designer_m).npv)
    else
        obj = mean(- compute_annualised_capex(1, s, des_m, designer_m) - compute_grid_cost(2, s, des_m) for s in 1:ns)
    end

    # Add the LPSP constraint for the heat
    isa(des_m.ld_H, Load) ? obj -= 1e32 * max(0., mean(LPSP(y, s, des_m).lpsp_H - 0.05 for s in 1:ns, y in 2:ny)) : nothing

    # Add the soc constraint for the seasonal storage
    isa(des_m.h2tank, H2Tank) ? obj -= 1e32 *  max(0., mean(des_m.h2tank.soc[1,y,s] - des_m.h2tank.soc[end,y,s] for s in 1:ns, y in 2:ny)) : nothing

    # Add the share constraint
    designer.options.share_constraint ? obj -= 1e32 * max(0., des.parameters.τ_share - mean(compute_share(y, s, des_m) for s in 1:ns, y in 2:ny)) : nothing

    return obj
end

### Offline
function initialize_designer!(des::DistributedEnergySystem, designer::Metaheuristic, ω::AbstractScenarios)
    # Preallocate and assigned values
    preallocate!(designer, des.parameters.ny, des.parameters.ns)

    # Scenario reduction from the optimization scenario pool
    if designer.options.scenario_reduction
        if designer.options.mode == "deterministic"
            ω_m = scenarios_reduction(ω, 1:des.parameters.nh, 1, 1)
            # Concatenation to simulate 2 years
            ω_m = concatenate(ω_m, ω_m, dims=2)
        elseif designer.options.mode == "twostage"
            ω_m = scenarios_reduction(ω, 1:des.parameters.nh, 2:2, 1:5)
            # Concatenation to simulate 2 years
            ω_m = concatenate(ω_m, ω_m, dims=2)
        elseif designer.options.mode == "npv"
            ω_m = scenarios_reduction(ω, 1:des.parameters.nh, 1:des.parameters.ny, 1:1)
        end
    end

    # Bounds
    lb, ub = set_bounds(des)

    # Optimize
    designer.results = Metaheuristics.optimize(lb, ub,
                                               designer.options.method,
                                               options = Metaheuristics.Options(iterations=designer.options.iterations, multithreads=true)
    ) do decisions
        fobj(decisions, des, designer, ω_m)
      end

    # Assign values
    designer.u.pv[1,:] .= designer.results.minimizer[1]
    designer.u.liion[1,:] .= designer.results.minimizer[2]
    designer.u.h2tank[1,:] .= designer.results.minimizer[3]
    designer.u.elyz[1,:] .= designer.results.minimizer[4]
    designer.u.fc[1,:] .= designer.results.minimizer[5]
    designer.u.tes[1,:] .= designer.results.minimizer[6]

    # Save history for online optimization
    designer.history = ω

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
function copy(des::DistributedEnergySystem, nh::Int64, ny::Int64, ns::Int64)
    des_copy = DistributedEnergySystem(ld_E = isa(des.ld_E, Load) ? Load() : nothing,
                                  ld_H = isa(des.ld_H, Load) ? Load() : nothing,
                                  pv = isa(des.pv, Source) ? Source() : nothing,
                                  liion = isa(des.liion, Liion) ? Liion() : nothing,
                                  tes = isa(des.tes, ThermalSto) ? ThermalSto() : nothing,
                                  h2tank = isa(des.h2tank, H2Tank) ? H2Tank() : nothing,
                                  elyz = isa(des.elyz, Electrolyzer) ? Electrolyzer() : nothing,
                                  fc = isa(des.fc, FuelCell) ? FuelCell() : nothing,
                                  heater = isa(des.heater, Heater) ? Heater() : nothing,
                                  grid = isa(des.grid, Grid) ? Grid() : nothing,
                                  parameters = Genesys.GlobalParameters(nh, ny, ns, τ_share = des.parameters.τ_share))

    return des_copy
end
