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

    Metaheuristic(; options = MetaheuristicOptions()) = new(options)
end

### Offline
function initialize_designer!(des::DistributedEnergySystem, designer::Metaheuristic, ω_optim::AbstractScenarios)

     # Scenario reduction from the optimization scenario pool
     ω_meta = scenarios_reduction(designer, ω_optim)

     # Define the objective function from global variables
     function fobj(decisions::Array{Float64,1})
         # Initialize DES
         des_meta = deepcopy(des)

         # Initialize controller
         controller_meta = initialize_controller!(des_meta, designer.options.controller, ω_meta)

         # Initialize with the dummy designer
         designer_meta = initialize_designer!(des_meta, DummyDesigner(), ω_meta)

         # Initialize with the decisions variables
         isa(des.pv, Source) ? designer_meta.u.pv[1] = decisions[1] : nothing
         isa(des.liion, Liion) ? designer_meta.u.liion[1] = decisions[2]  : nothing
         isa(des.h2tank, H2Tank) ? designer_meta.u.h2tank[1] = decisions[3]  : nothing
         isa(des.elyz, Electrolyzer) ? designer_meta.u.elyz[1] = decisions[4]  : nothing
         isa(des.fc, FuelCell) ? designer_meta.u.fc[1] = decisions[5]  : nothing
         isa(des.tes, ThermalSto) ? designer_meta.u.tes[1] = decisions[6]  : nothing

         # Simulate
         simulate!(des_meta, controller_meta, designer_meta, ω_meta)

         # Compute cost indicators
         costs_meta = compute_economics(des_meta, designer_meta)

         # Compute tech indicators
         tech_meta = compute_tech_indicators(des_meta)

         # Objective - algorithm find the maximum
        designer.option.share_constraint ? obj = costs_meta.npv[1] - 1e32 * max(0., des.parameters.τ_share - minimum(tech_meta.τ_share[2:end,:])) : obj = costs_meta.npv[1]

         return obj
     end

     # Bounds
     lb, ub = Float64[], Float64[]
     if isa(des.pv, Source)
         push!(lb, 0.) ; push!(ub, 1000.)
     end
     if isa(des.liion, Liion)
         push!(lb, 0.) ; push!(ub, 1000.)
     end
     if isa(des.h2tank, H2Tank)
         push!(lb, 0.) ; push!(ub, 50000.)
     end
     if isa(des.elyz, Electrolyzer)
         push!(lb, 0.) ; push!(ub, 50.)
     end
     if isa(des.fc, FuelCell)
         push!(lb, 0.) ; push!(ub, 50.)
     end
     if isa(des.tes, ThermalSto)
         push!(lb, 0.) ; push!(ub, 1000.)
     end

     # Optimize
     designer.results = Metaheuristics.optimize(fobj, lb, ub,
                                                designer.options.method,
                                                options = Metaheuristics.Options(iterations=designer.options.iterations, multithreads=true))

     # Preallocate and assigned values
     preallocate!(designer, des.parameters.ny, des.parameters.ns)
     isa(des.pv, Source) ? designer.u.pv[1,:] .= designer.results.minimizer[1] : nothing
     isa(des.liion, Liion) ? designer.u.liion[1,:] .= designer.results.minimizer[2] : nothing
     isa(des.h2tank, H2Tank) ? designer.u.h2tank[1,:] .= designer.results.minimizer[3] : nothing
     isa(des.elyz, Electrolyzer) ? designer.u.elyz[1,:] .= designer.results.minimizer[4] : nothing
     isa(des.fc, FuelCell) ? designer.u.fc[1,:] .= designer.results.minimizer[5] : nothing
     isa(des.tes, ThermalSto) ? designer.u.tes[1,:] .= designer.results.minimizer[5] : nothing

     return designer
end

### Online
function compute_investment_decisions!(y::Int64, s::Int64, des::DistributedEnergySystem, designer::Metaheuristic)
    ϵ = 0.1

    if designer.options.reopt
        # Do we need to reoptimize ?
        (isa(des.liion, Liion) && des.liion.soh[end,y,s] < ϵ) || (isa(des.elyz, Electrolyzer) && des.elyz.soh[end,y,s] < ϵ) || (isa(des.fc, FuelCell) && des.fc.soh[end,y,s] < ϵ) ? nothing : return
        println("Re-optimization not yet implemented...")
    else
        isa(des.liion, Liion) && des.liion.soh[end,y,s] < ϵ ? designer.u.liion[y,s] = des.liion.Erated[y,s] : nothing
        isa(des.elyz, Electrolyzer) && des.elyz.soh[end,y,s] < ϵ ? designer.u.elyz[y,s] = des.elyz.powerMax[y,s] : nothing
        isa(des.fc, FuelCell) && des.fc.soh[end,y,s] < ϵ ? designer.u.fc[y,s] = des.fc.powerMax[y,s] : nothing
    end
end
