#=
    Anticipative multi-stage designer with anticipative controller
=#

mutable struct AnticipativeMultiStageDesigner <: AbstractMultiStageDesigner
    u::NamedTuple
    model::JuMP.Model
    parameters::Dict{String, Any}
    AnticipativeMultiStageDesigner() = new()
end

#### Models ####
# Simple
function multistage_milp_model(ld::Load, pv::Source, liion::Liion,
    controller::AnticipativeController, designer::AnticipativeMultiStageDesigner,
     grid::Grid, ω_optim::Scenarios, parameters::NamedTuple)

    # Parameters
    M = 1000 # big-M value
    γ = 1. ./ (1. + parameters.τ) .^ (1:parameters.Y) # discount factor

    # Sets
    nh = size(ω_optim.values.ld_E,1) # Number of hours
    ny = size(ω_optim.values.ld_E,2) # Number of years

    # Model definition
    m = Model(CPLEX.Optimizer)

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
    [h in 1:nh, y in 1:ny], p_liion_dch[h,y] <= liion.α_p_dch * E_liion[y]
    [h in 1:nh, y in 1:ny], p_liion_ch[h,y] >= -liion.α_p_ch * E_liion[y]
    [h in 1:nh, y in 2:ny], p_g_in[h,y] <= (1. - grid.τ_power) * maximum(ω_optim.values.ld_E)
    # Operation state bounds
    [h in 1:nh+1, y in 1:ny], soc_liion[h,y] <= liion.α_soc_max * E_liion[y]
    [h in 1:nh+1, y in 1:ny], soc_liion[h,y] >= liion.α_soc_min * E_liion[y]
    [h in 1:nh+1, y in 1:ny], soh_liion[h,y] <= E_liion[y]
    # Operation state dynamics
    [h in 1:nh, y in 1:ny], soc_liion[h+1,y] == soc_liion[h,y] * (1. - liion.η_self * parameters.Δh) - (p_liion_ch[h,y] * liion.η_ch + p_liion_dch[h,y] / liion.η_dch) * parameters.Δh
    [h in 1:nh, y in 1:ny], soh_liion[h+1,y] == soh_liion[h,y] - (p_liion_dch[h,y] - p_liion_ch[h,y]) * parameters.Δh / (2. * liion.dod * liion.nCycle)
    # Power balance
    [h in 1:nh, y in 1:ny], ω_optim.values.ld_E[h,y] - pMax_pv[y] *  ω_optim.values.pv_E[h,y] <= p_g_out[h,y] + p_g_in[h,y] + p_liion_ch[h,y] + p_liion_dch[h,y]
    # Self-sufficiency constraint
    [y in 2:ny], sum(p_g_in[h,y] for h in 1:nh) <= (1. - grid.τ_energy) * sum(ω_optim.values.ld_E[h,y] for h in 1:nh)
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
    [y in 1:ny-1], soc_liion[1,y+1] - liion.α_soc_max * r_liion[y] <= M * (1 - δ_inv_liion[y])
    [y in 1:ny-1], liion.α_soc_max * r_liion[y] - soc_liion[1,y+1] <= M * (1 - δ_inv_liion[y])
    # SoH continuity
    [y in 1:ny-1], soh_liion[1,y+1] - soh_liion[nh+1,y] <= M * δ_inv_liion[y]
    [y in 1:ny-1], soh_liion[nh+1,y] - soh_liion[1,y+1] <= M * δ_inv_liion[y]
    [y in 1:ny-1], soh_liion[1,y+1] - r_liion[y] <= M * (1 - δ_inv_liion[y])
    [y in 1:ny-1], r_liion[y] - soh_liion[1,y+1] <= M * (1 - δ_inv_liion[y])
    # Initial and final conditions
    E_liion[1] == liion.Erated[1]
    pMax_pv[1] == pv.powerMax[1]
    soh_liion[1,1] == liion.soh[1,1] * liion.Erated[1]
    soc_liion[1,1] == liion.soc[1,1] * liion.Erated[1]
    [y in 1:ny], soc_liion[1,y] == soc_liion[nh+1,y]
    end)

    # Objective
    @objective(m, Min, sum( γ[y] * ( ω_optim.values.C_liion[y] * r_liion[y] + ω_optim.values.C_pv[y] * r_pv[y] +
    sum( (p_g_in[h, y] * ω_optim.values.C_grid_in[h,y] + p_g_out[h, y] * ω_optim.values.C_grid_out[h,y]) * parameters.Δh for h=1:nh) ) for y=1:ny))

    return m
end

#### Offline functions ####
# Simple
function offline_optimization(ld::Load, pv::Source, liion::Liion,
    controller::AnticipativeController, designer::AnticipativeMultiStageDesigner,
     grid::Grid, ω_optim::Scenarios, parameters::NamedTuple)

     # Scenario reduction from the optimization scenario pool
     ω_anticipative = scenarios_reduction(designer, ω_optim)

     # Initialize model
     designer.model = controller.model = multistage_milp_model(ld, pv, liion, controller, designer, grid, ω_anticipative, parameters)

     # Compute both investment and operation decisions
     optimize!(designer.model)

     # Formatting variables to simulation
     # Operation decisions
     controller.u = (
     u_liion =  value.(controller.model[:p_liion_ch] + controller.model[:p_liion_dch]),
     )
     # Investment decisions
     designer.u =(
     u_liion = value.(designer.model[:r_liion]),
     u_pv = value.(designer.model[:r_pv]),
     )
end

#### Online functions ####
# Simple
function compute_investment_decisions(y::Int64, s::Int64, ld::Load, pv::Source,
    liion::Liion, grid::Grid, designer::AnticipativeMultiStageDesigner, ω_optim::Scenarios, parameters::NamedTuple)
    return nothing
end
