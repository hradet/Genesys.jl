# Anticipative controller

mutable struct AnticipativeOptions
    solver

    AnticipativeOptions(; solver = CPLEX) = new(solver)
end

mutable struct Anticipative <: AbstractController
    options::AnticipativeOptions
    generations::Vector{Float64}
    storages::Vector{Float64}
    converters::Vector{Float64}
    decisions::NamedTuple

    Anticipative(; options = AnticipativeOptions(),
                   generations = [0.],
                   storages = [0.],
                   converters = [0.]) =
                   new(options, generations, storages, converters)
end

### Models
function build_model(mg::Microgrid, controller::Anticipative, ω::Scenarios)
    # Initialize
    m = Model(controller.options.solver.Optimizer)
    set_optimizer_attribute(m,"CPX_PARAM_SCRIND", 0)
    # Add investment variables
    add_investment_decisions!(m, mg, ω)
    # Fix their values
    fix_investment_decisions!(m, mg, controller)
    # Add decision variables
    add_operation_decisions!(m, mg, ω)
    # Add technical constraints
    add_technical_constraints!(m, mg, ω)
    # Add power balance constraints
    add_power_balances!(m, mg, ω)
    # Objective
    opex = compute_opex(m, mg, ω)
    @objective(m, Min, opex[1])
    return m
end

### Offline
function initialize_controller!(mg::Microgrid, controller::Anticipative, ω::Scenarios)
    # Preallocate
    preallocate!(mg, controller)

    for y in 2:mg.parameters.ny, s in 1:mg.parameters.ns
        # Scenario reduction
        ω_reduced, _ = reduce(ManualReducer(h = 1:mg.parameters.nh, y = y:y, s = s:s), ω)
        # Build model
        model = build_model(mg, controller, ω_reduced)
        # Optimize
        optimize!(model)
        # Assign controller values
        for k in 1:length(mg.storages)
            controller.decisions.storages[k][:,y,s] .= value.(model[:p_dch][:,1,k] .- model[:p_ch][:,1,k])
        end
        for (k,a) in enumerate(mg.converters)
            if a isa Heater
                controller.decisions.converters[k][:,y,s] .= .- value.(model[:p_c][:,1,k])
            elseif a isa Electrolyzer
                controller.decisions.converters[k][:,y,s] .= .- value.(model[:p_c][:,1,k])
            elseif a isa FuelCell
                controller.decisions.converters[k][:,y,s] .= value.(model[:p_c][:,1,k])
            end
        end
     end

     return controller
end

### Online
function compute_operation_decisions!(h::Int64, y::Int64, s::Int64, mg::Microgrid, controller::Anticipative)
    return controller
end
