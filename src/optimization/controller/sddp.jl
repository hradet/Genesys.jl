#=
    Stochastic dual dynamic programming controller
=#

mutable struct SDDPCOptions
    solver
    reducer::AbstractScenariosReducer
    risk::SDDP.AbstractRiskMeasure
    parallel::SDDP.AbstractParallelScheme
    iterations::Int64
    seasonal_targets::Union{Array{Float64,2}, Nothing}
    read_cuts::Union{String, Nothing}
    write_cuts::Union{String, Nothing}

    SDDPCOptions(; solver = CPLEX,
                  reducer = FeatureBasedReducer(),
                  risk = SDDP.Expectation(),
                  parallel = SDDP.Serial(),
                  iterations = 100,
                  seasonal_targets = nothing,
                  read_cuts = nothing,
                  write_cuts = nothing) = new(solver, reducer, risk, parallel, iterations, seasonal_targets, read_cuts, write_cuts)
end

# SDDP controller definition
mutable struct SDDPC <: AbstractController
    options::SDDPCOptions
    pv::Float64
    liion::Float64
    tes::Float64
    h2tank::Float64
    elyz::Float64
    fc::Float64
    u::NamedTuple
    policy::Array{SDDP.DecisionRule{Int64}, 1}
    model::SDDP.PolicyGraph{Int64}
    SDDPC(; options = SDDPCOptions(),
                   pv = 0.,
                   liion = 0.,
                   tes = 0.,
                   h2tank = 0.,
                   elyz = 0.,
                   fc = 0.) =
                   new(options, pv, liion, tes, h2tank, elyz, fc)
end

### Model
function build_model(des::DistributedEnergySystem, controller::SDDPC, ω::Scenarios, probabilities::Vector{Float64})
    # Parameters
    λ = 1e3
    nh = size(ω.ld_E.power, 1) # Stages
    ns = size(ω.ld_E.power, 3) # Scenarios

    # Initialize
    isa(des.liion, Liion) ? liion = des.liion : liion = Liion()
    isa(des.tes, ThermalSto) ? tes = des.tes : tes = ThermalSto()
    isa(des.h2tank, H2Tank) ? h2tank = des.h2tank : h2tank = H2Tank()
    isa(des.elyz, Electrolyzer) ? elyz = des.elyz : elyz = Electrolyzer()
    isa(des.fc, FuelCell) ? fc = des.fc : fc = FuelCell()
    isa(des.heater, Heater) ? heater = des.heater : heater = Heater()
    isa(des.pv, Source) ? pv = des.pv : pv = Source()
    isa(des.grid, Grid) ? grid = des.grid : grid = Grid()

    # Model
    model = SDDP.LinearPolicyGraph(
        stages = nh, sense = :Min, lower_bound = -1000., optimizer = controller.options.solver.Optimizer
    ) do sp, h
        # Define the state variable.
        @variable(sp, liion.α_soc_min * controller.liion <= soc_liion <= liion.α_soc_max * controller.liion, SDDP.State, initial_value = liion.soc_ini * controller.liion)
        @variable(sp, tes.α_soc_min * controller.tes <= soc_tes <= tes.α_soc_max * controller.tes, SDDP.State, initial_value = tes.soc_ini * controller.tes)
        @variable(sp, h2tank.α_soc_min * controller.h2tank <= soc_h2tank <= h2tank.α_soc_max * controller.h2tank, SDDP.State, initial_value = h2tank.soc_ini * controller.h2tank)
        # Define the control variables.
        @variables(sp, begin
            p_liion_ch      <= 0.
            p_liion_dch     >= 0.
            p_tes_ch        <= 0.
            p_tes_dch       >= 0.
            p_h2tank_ch     <= 0.
            p_h2tank_dch    >= 0.
            p_elyz_E        <= 0.
            p_fc_E          >= 0.
            p_heater_E      <= 0.
            p_g_out         <= 0.
            p_g_in          >= 0.
            # Variables to be parametrized
            ld_E
            ld_H
            pv
            # Deterministic target penalty - aux variable for max linearization
            target          >= 0.
        end)
        # Define the constraints
        @constraints(sp, begin
            # Power bounds
            p_liion_dch  <= liion.α_p_dch * controller.liion
            p_liion_ch   >= -liion.α_p_ch * controller.liion
            p_tes_dch    <= tes.α_p_dch * controller.tes
            p_tes_ch     >= -tes.α_p_ch * controller.tes
            p_h2tank_dch <= h2tank.α_p_dch * controller.h2tank
            p_h2tank_ch  >= -h2tank.α_p_ch * controller.h2tank
            p_elyz_E     >= -controller.elyz
            p_fc_E       <= controller.fc
            p_heater_E   >= -heater.powerMax_ini
            p_g_in       <= grid.powerMax
            p_g_out      >= -grid.powerMax
            # Dynamics
            soc_liion.out  == soc_liion.in * (1. - liion.η_self * des.parameters.Δh) - (p_liion_ch * liion.η_ch + p_liion_dch / liion.η_dch) * des.parameters.Δh
            soc_tes.out    == soc_tes.in * (1. - tes.η_self * des.parameters.Δh) - (p_tes_ch * tes.η_ch + p_tes_dch / tes.η_dch) * des.parameters.Δh
            soc_h2tank.out == soc_h2tank.in * (1. - h2tank.η_self * des.parameters.Δh) - (p_h2tank_ch * h2tank.η_ch + p_h2tank_dch / h2tank.η_dch) * des.parameters.Δh
            # Power balances
            ld_E <= controller.pv * pv + p_liion_ch + p_liion_dch + p_elyz_E + p_fc_E + p_heater_E + p_g_in + p_g_out
            ld_H <= p_tes_ch  + p_tes_dch - elyz.η_E_H * p_elyz_E + fc.η_H2_H / fc.η_H2_E * p_fc_E - heater.η_E_H * p_heater_E
            0.   == p_h2tank_ch + p_h2tank_dch - elyz.η_E_H2 * p_elyz_E - p_fc_E / fc.η_H2_E
        end)
        # Deterministic targets for the seasonal storage
        if !isa(controller.options.seasonal_targets, Nothing)
            @constraints(sp, begin
            target >= controller.options.seasonal_targets[h] - soc_h2tank.out
            end)
        end
        # Parameterize the subproblem.
        Ω = [(ld_E = ω.ld_E.power[h,1,s], ld_H = ω.ld_H.power[h,1,s], pv = ω.pv.power[h,1,s], cost_in = ω.grid.cost_in[h,1,s], cost_out = ω.grid.cost_out[h,1,s]) for s in 1:ns]
        SDDP.parameterize(sp, Ω, probabilities) do ω
            JuMP.fix(ld_E, ω.ld_E)
            JuMP.fix(ld_H, ω.ld_H)
            JuMP.fix(pv, ω.pv)
            @stageobjective(sp, p_g_in * ω.cost_in + p_g_out * ω.cost_out + λ * target)
        end
    end

    return model
end

### Offline
function initialize_controller!(des::DistributedEnergySystem, controller::SDDPC, ω::Scenarios)
    # Scenario reduction to train the model
    println("Starting scenario reduction...")
    ω_reduced, probabilities = reduce(controller.options.reducer, ω)

    # Build model
    println("Building the model...")
    controller.model = build_model(des, controller, ω_reduced, probabilities)

    # Train the model or read cuts from file
    if isa(controller.options.read_cuts, Nothing)
        # Train policy
        println("Starting offline training...")
        SDDP.train(controller.model,
                   risk_measure = controller.options.risk,
                   parallel_scheme = controller.options.parallel,
                   iteration_limit = controller.options.iterations,
                   print_level = 0)
        # Save cuts to file
        if !isa(controller.options.write_cuts, Nothing)
            println("Writting cuts to file...")
            SDDP.write_cuts_to_file(controller.model, controller.options.write_cuts)
        end
    else
        println("Reading cuts from file...")
        SDDP.read_cuts_from_file(controller.model, controller.options.read_cuts)
    end

    # Store policy for each stage
    controller.policy = [SDDP.DecisionRule(controller.model, node = h) for h in 1:des.parameters.nh]

    # Preallocation
    preallocate!(controller, des.parameters.nh, des.parameters.ny, des.parameters.ns)

    return controller
end

### Online
function compute_operation_decisions!(h::Int64, y::Int64, s::Int64, des::DistributedEnergySystem, controller::SDDPC)
    # Evaluate the policy at each stage
    _, _, controls = SDDP.evaluate(controller.policy[h],
                      incoming_state = Dict(:soc_liion => des.liion.soc[h,y,s] * des.liion.Erated[y,s], :soc_tes => des.tes.soc[h,y,s] * des.tes.Erated[y,s], :soc_h2tank => des.h2tank.soc[h,y,s] * des.h2tank.Erated[y,s]),
                      noise = (ld_E = des.ld_E.power[h,y,s], ld_H = des.ld_H.power[h,y,s], pv = des.pv.power_E[h,y,s] / des.pv.powerMax[y,s], cost_in = des.grid.cost_in[h,y,s], cost_out = des.grid.cost_out[h,y,s]),
                      controls_to_record = [:p_liion_ch, :p_liion_dch, :p_tes_ch, :p_tes_dch, :p_h2tank_ch, :p_h2tank_dch, :p_elyz_E, :p_fc_E, :p_heater_E])

    # Operation decision
    controller.u.liion[h,y,s] = controls[:p_liion_ch] + controls[:p_liion_dch]
    controller.u.tes[h,y,s] = controls[:p_tes_ch] + controls[:p_tes_dch]
    controller.u.h2tank[h,y,s] = controls[:p_h2tank_ch] + controls[:p_h2tank_dch]
    controller.u.elyz[h,y,s] = controls[:p_elyz_E]
    controller.u.fc[h,y,s] = controls[:p_fc_E]
    controller.u.heater[h,y,s] = controls[:p_heater_E]

    return controller
end
