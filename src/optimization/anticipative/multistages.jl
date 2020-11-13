#=
    Anticipative multi-stage designer with anticipative controller
=#

mutable struct AnticipativeMultiStagesOptions
    solver
    scenario_reduction::String
    s::Int64

    AnticipativeMultiStagesOptions(; solver = CPLEX, scenario_reduction = "manual", s = 1) = new(solver, scenario_reduction, s)
end

mutable struct AnticipativeMultiStages <: AbstractDesigner
    options::AnticipativeMultiStagesOptions
    u::NamedTuple
    model::JuMP.Model
    AnticipativeMultiStages(; options = AnticipativeMultiStagesOptions()) = new(options)
end

### Models
function build_model(des::DistributedEnergySystem, designer::AnticipativeMultiStages, ω::AbstractScenarios)

    # Parameters
    M = 1000 # big-M value
    γ = 1. ./ (1. + des.parameters.τ) .^ (1:des.parameters.ny) # discount factor

    # Sets
    nh = des.parameters.nh # Number of hours
    ny = des.parameters.ny # Number of years

    # Model definition
    m = Model(designer.options.solver.Optimizer)

    if isa(des.liion, Liion)
        # Variables
        @variables(m, begin
        # Operation decision variables
        p_liion_ch[1:nh, 1:ny] <= 0.
        p_liion_dch[1:nh, 1:ny] >= 0.
        p_g_out[1:nh, 1:ny] <= 0.
        p_g_in[1:nh, 1:ny] >= 0.
        # Investment decision variables
        r_liion[1:ny] >= 0.
        δ_inv_liion[1:ny], (Bin)
        r_pv[1:ny] >= 0.
        δ_inv_pv[1:ny], (Bin)
        # Operation state variables
        soc_liion[1:nh+1, 1:ny]
        soh_liion[1:nh+1, 1:ny] >= 0.
        # Investment state variables
        0 <= E_liion[1:ny+1] <= 1000
        0 <= pMax_pv[1:ny+1] <= 1000
        end)
        # Constraints
        @constraints(m, begin
        # Power bounds
        [h in 1:nh, y in 1:ny], p_liion_dch[h,y] <= des.liion.α_p_dch * E_liion[y]
        [h in 1:nh, y in 1:ny], p_liion_ch[h,y] >= -des.liion.α_p_ch * E_liion[y]
        [h in 1:nh, y in 2:ny], p_g_in[h,y] <= des.grid.powerMax
        # Operation state bounds
        [h in 1:nh+1, y in 1:ny], soc_liion[h,y] <= des.liion.α_soc_max * E_liion[y]
        [h in 1:nh+1, y in 1:ny], soc_liion[h,y] >= des.liion.α_soc_min * E_liion[y]
        [h in 1:nh+1, y in 1:ny], soh_liion[h,y] <= E_liion[y]
        # Operation state dynamics
        [h in 1:nh, y in 1:ny], soc_liion[h+1,y] == soc_liion[h,y] * (1. - des.liion.η_self * des.parameters.Δh) - (p_liion_ch[h,y] * des.liion.η_ch + p_liion_dch[h,y] / des.liion.η_dch) * des.parameters.Δh
        [h in 1:nh, y in 1:ny], soh_liion[h+1,y] == soh_liion[h,y] - (p_liion_dch[h,y] - p_liion_ch[h,y]) * des.parameters.Δh / (2. * (des.liion.α_soc_max - des.liion.α_soc_min) * des.liion.nCycle)
        # Power balance
        [h in 1:nh, y in 1:ny], ω.ld_E.power[h,y] - pMax_pv[y] *  ω.pv.power[h,y] <= p_g_out[h,y] + p_g_in[h,y] + p_liion_ch[h,y] + p_liion_dch[h,y]
        # Self-sufficiency constraint
        [y in 2:ny], sum(p_g_in[h,y] for h in 1:nh) <= (1. - des.parameters.τ_share) * sum(ω.ld_E.power[h,y] for h in 1:nh)
        # Investment bounds
        [y in 1:ny], r_liion[y] <= 1000. * δ_inv_liion[y]
        [y in 1:ny], r_pv[y] <= 1000. * δ_inv_pv[y]
        # Investment state dynamics
        # if δ_inv_liion = 0 (r_liion = 0) then E_liion[y+1] =  E_liion[y]
        [y in 1:ny], E_liion[y+1] - E_liion[y] <= M * δ_inv_liion[y]
        [y in 1:ny], E_liion[y] - E_liion[y+1] <= M * δ_inv_liion[y]
        # else E_liion[y+1] = r_liion[y] end
        [y in 1:ny], E_liion[y+1] - r_liion[y] <= M * (1 - δ_inv_liion[y])
        [y in 1:ny], r_liion[y] - E_liion[y+1] <= M * (1 - δ_inv_liion[y])
        # PV
        [y in 1:ny], pMax_pv[y+1] - pMax_pv[y] <= M * δ_inv_pv[y]
        [y in 1:ny], pMax_pv[y] - pMax_pv[y+1] <= M * δ_inv_pv[y]
        [y in 1:ny], pMax_pv[y+1] - r_pv[y] <= M * (1 - δ_inv_pv[y])
        [y in 1:ny], r_pv[y] - pMax_pv[y+1] <= M * (1 - δ_inv_pv[y])
        # SoC continuity
        [y in 1:ny-1], soc_liion[1,y+1] - soc_liion[nh+1,y] <= M * δ_inv_liion[y]
        [y in 1:ny-1], soc_liion[nh+1,y] - soc_liion[1,y+1] <= M * δ_inv_liion[y]
        [y in 1:ny-1], soc_liion[1,y+1] - des.liion.α_soc_max * r_liion[y] <= M * (1 - δ_inv_liion[y])
        [y in 1:ny-1], des.liion.α_soc_max * r_liion[y] - soc_liion[1,y+1] <= M * (1 - δ_inv_liion[y])
        # SoH continuity
        [y in 1:ny-1], soh_liion[1,y+1] - soh_liion[nh+1,y] <= M * δ_inv_liion[y]
        [y in 1:ny-1], soh_liion[nh+1,y] - soh_liion[1,y+1] <= M * δ_inv_liion[y]
        [y in 1:ny-1], soh_liion[1,y+1] - r_liion[y] <= M * (1 - δ_inv_liion[y])
        [y in 1:ny-1], r_liion[y] - soh_liion[1,y+1] <= M * (1 - δ_inv_liion[y])
        # Initial and final conditions
        E_liion[1] == des.liion.Erated[1,1]
        pMax_pv[1] == des.pv.powerMax[1,1]
        soh_liion[1,1] == des.liion.soh[1,1,1] * des.liion.Erated[1,1]
        soc_liion[1,1] == des.liion.soc[1,1,1] * des.liion.Erated[1,1]
        [y in 1:ny], soc_liion[1,y] == soc_liion[nh+1,y]
        end)
    end

    # Objective
    @objective(m, Min, sum( γ[y] * ( ω.liion.cost[y] * r_liion[y] + ω.pv.cost[y] * r_pv[y] +
    sum( (p_g_in[h, y] * ω.grid.cost_in[h,y] + p_g_out[h, y] * ω.grid.cost_out[h,y]) * des.parameters.Δh for h=1:nh) ) for y=1:ny))

    return m
end

### Offline
function offline_optimization!(des::DistributedEnergySystem, designer::AnticipativeMultiStages, ω::AbstractScenarios)

     # Scenario reduction from the optimization scenario pool
     ω_anticipative = scenarios_reduction(designer, ω)

     # Build model
     designer.model = build_model(des, designer, ω_anticipative)

     # Compute both investment and operation decisions
     optimize!(designer.model)

     # Preallocate and assigned values to the controller
     controller = Anticipative()
     preallocate!(controller, des.parameters.nh, des.parameters.ny, des.parameters.ns)
     isa(des.liion, Liion) ? controller.u.liion .= value.(designer.model[:p_liion_ch] + designer.model[:p_liion_dch]) : nothing

     # Preallocate and assigned values to the designer
     preallocate!(designer, des.parameters.ny, des.parameters.ns)
     isa(des.liion, Liion) ? designer.u.liion .= value.(designer.model[:r_liion]) : nothing
     isa(des.pv, Source) ? designer.u.pv .= value.(designer.model[:r_pv]) : nothing

     return controller, designer
end

### Online
function compute_investment_decisions!(y::Int64, s::Int64, des::DistributedEnergySystem, designer::AnticipativeMultiStages)
    return designer
end
