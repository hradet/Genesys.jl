#=
    Designer based on the equivalent annual cost (EAC) with multiple scenarios
=#

mutable struct MILPOptions
  solver::Module
  reducer::AbstractScenariosReducer
  objective_risk::AbstractRiskMeasure
  share_risk::AbstractRiskMeasure
  reopt::Bool
  read_reduction::Union{String, Nothing}
  write_reduction::Union{String, Nothing}

  MILPOptions(; solver = Cbc,
                reducer = FeatureBasedReducer(),
                objective_risk = Expectation(),
                share_risk = Expectation(),
                reopt = false,
                read_reduction = nothing,
                write_reduction = nothing) =
                new(solver, reducer, objective_risk, share_risk, reopt, read_reduction, write_reduction)
end

mutable struct MILP <: AbstractDesigner
    options::MILPOptions
    decisions::NamedTuple
    model::JuMP.Model
    history::AbstractScenarios

    MILP(; options = MILPOptions()) = new(options)
end

### Models
function build_model(mg::Microgrid, designer::MILP, ω::Scenarios, probabilities::Vector{Float64})
    # Sets
    nh, ns = size(ω.demands[1].power, 1), size(ω.demands[1].power, 3)
    # Initialize
    m = Model(designer.options.solver.Optimizer)
    # Add desgin decision variables
    add_investment_decisions!(m, mg.generations)
    add_investment_decisions!(m, mg.storages)
    add_investment_decisions!(m, mg.converters)
    # Add operation decision variables
    add_operation_decisions!(m, mg.storages, nh, ns)
    add_operation_decisions!(m, mg.converters, nh, ns)
    add_operation_decisions!(m, mg.grids, nh, ns)
    # Add design constraints
    add_investment_constraints!(m, mg.generations)
    add_investment_constraints!(m, mg.storages)
    add_investment_constraints!(m, mg.converters)
    # Add technical constraints
    add_technical_constraints!(m, mg.storages, mg.parameters.Δh, nh, ns)
    add_technical_constraints!(m, mg.converters, nh, ns)
    add_technical_constraints!(m, mg.grids, nh, ns)
    # Add periodicity constraint
    add_periodicity_constraints!(m, mg.storages, ns)
    # Add power balance constraints
    add_power_balance!(m, mg, ω, Electricity, nh, ns)
    add_power_balance!(m, mg, ω, Heat, nh, ns)
    add_power_balance!(m, mg, ω, Hydrogen, nh, ns)
    # Renewable share constraint
    add_renewable_share!(m, mg, ω, probabilities, designer.options.share_risk, nh, ns)
    # Objective
    add_design_objective!(m, mg, ω, probabilities, designer.options.objective_risk, nh, ns)
    return m
end

### Offline
function initialize_designer!(mg::Microgrid, designer::MILP, ω::Scenarios)
    # Preallocate
    preallocate!(mg, designer)

    # Scenario reduction from the optimization scenario pool
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

    # Initialize model
    println("Building the model...")
    designer.model = build_model(mg, designer, ω_reduced, probabilities)

    # Compute investment decisions for the first year
    println("Starting optimization...")
    optimize!(designer.model)

    # Assign values
    for k in 1:length(mg.generations)
        designer.decisions.generations[k][1,:] .= value(designer.model[:r_g][k])
    end
    for k in 1:length(mg.storages)
        designer.decisions.storages[k][1,:] .= value(designer.model[:r_sto][k])
    end
    for k in 1:length(mg.converters)
        designer.decisions.converters[k][1,:] .= value(designer.model[:r_c][k])
    end

    # Save history
    designer.history = ω_reduced

     return designer
end

### Online
function compute_investment_decisions!(y::Int64, s::Int64, mg::Microgrid, designer::MILP)
    # TODO
end

### Utils
beta(risk::WorstCase) = 1. - 1e-6
beta(risk::Expectation) = 0.
beta(risk::CVaR) = risk.β
Γ(τ::Float64, lifetime::Union{Float64, Int64}) = τ * (τ + 1.) ^ lifetime / ((τ + 1.) ^ lifetime - 1.)
