#=
    Anticipative multi-stage designer with anticipative controller
=#

#                                   Model definition
#_______________________________________________________________________________
function multistage_milp_model(ld::Load, pv::Source, liion::Liion,
    controller::AnticipativeController, designer::AnticipativeMultiStageDesigner,
     grid::Grid, ω_tb::ClusteredScenarios, parameters::NamedTuple)

    # Parameters
    M = 1000 # big-M value
    γ = 1. ./ (1. + parameters.τ) .^ (1:parameters.Y) # discount factor

    # Sets
    ntb = size(ω_tb.clusters.ld_E,2) # Number of time block
    nh = length(1:parameters.Δh:controller.horizon) # Number of hours

    # Model definition
    m = Model(CPLEX.Optimizer)

    # Variables
    # Operation control variables
    @variables(m, begin
    # Liion
    p_liion_ch[1:nh, 1:ntb] <= 0.
    p_liion_dch[1:nh, 1:ntb] >= 0.
    # Recourse
    p_g_out[1:nh, 1:ntb] <= 0.
    p_g_in[1:nh, 1:ntb] >= 0.
    end)
    # Investment control variables
    @variables(m, begin
    r_liion[1:ntb] >= 0.
    δ_inv_liion[1:ntb], (Bin)
    r_pv[1:ntb] >= 0.
    δ_inv_pv[1:ntb], (Bin)
    end)
    # Operation state variables
    @variables(m, begin
    soc_liion[1:nh+1, 1:ntb]
    soh_liion[1:nh+1, 1:ntb] >= 0.
    end)
    # Investment state variables
    @variables(m, begin
    E_liion[1:ntb+1] >= 0.
    pMax_pv[1:ntb+1] >= 0.
    end)

    # Operation constraints bounds
    @constraints(m, begin
    # Controls
    # Liion
    [h in 1:nh, tb in 1:ntb], p_liion_dch[h,tb] <= liion.α_p_dch * E_liion[tb]
    [h in 1:nh, tb in 1:ntb], p_liion_ch[h,tb] >= -liion.α_p_ch * E_liion[tb]

    # State
    # Liion
    [h in 1:nh+1, tb in 1:ntb], soc_liion[h,tb] <= liion.α_soc_max * E_liion[tb]
    [h in 1:nh+1, tb in 1:ntb], soc_liion[h,tb] >= liion.α_soc_min * E_liion[tb]
    [h in 1:nh+1, tb in 1:ntb], soh_liion[h,tb] <= E_liion[tb]
    end)
    # Investment constraints bounds
    @constraints(m, begin
    # Controls
    [tb in 1:ntb], r_liion[tb] <= 1000. * δ_inv_liion[tb]
    [tb in 1:ntb], r_pv[tb] <= 1000. * δ_inv_pv[tb]
    end)
    # Operation constraints dynamics and recourse
    @constraints(m, begin
    # Dynamics
    # Liion
    [h in 1:nh, tb in 1:ntb], soc_liion[h+1,tb] == soc_liion[h,tb] * (1. - liion.η_self * parameters.Δh) - (p_liion_ch[h,tb] * liion.η_ch + p_liion_dch[h,tb] / liion.η_dch) * parameters.Δh
    [h in 1:nh, tb in 1:ntb], soh_liion[h+1,tb] == soh_liion[h,tb] -  ω_tb.nby[tb] * (p_liion_dch[h,tb] - p_liion_ch[h,tb]) * parameters.Δh / (2. * liion.dod * liion.nCycle)

    # Recourse
    [h in 1:nh, tb in 1:ntb], ω_tb.clusters.ld_E[h,tb] - pMax_pv[tb] *  ω_tb.clusters.pv_E[h,tb] <= p_g_out[h,tb] + p_g_in[h,tb] + p_liion_ch[h,tb] + p_liion_dch[h,tb]
    end)
    # Investment constraints dynamics
    @constraints(m, begin
    # Dynamics
    # E_liion
    # if δ_inv_liion = 0 (r_liion = 0) then E_liion[tb+1] =  E_liion[tb]
    [tb in 1:ntb], E_liion[tb+1] - E_liion[tb] <= M * δ_inv_liion[tb]
    [tb in 1:ntb], E_liion[tb] - E_liion[tb+1] <= M * δ_inv_liion[tb]
    # else E_liion[tb+1] = r_liion[tb] end
    [tb in 1:ntb], E_liion[tb+1] - r_liion[tb] <= M * (1 - δ_inv_liion[tb])
    [tb in 1:ntb], r_liion[tb] - E_liion[tb+1] <= M * (1 - δ_inv_liion[tb])
    # SOC_liion
    [tb in 1:ntb-1], soc_liion[1,tb+1] - soc_liion[nh+1,tb] <= M * δ_inv_liion[tb]
    [tb in 1:ntb-1], soc_liion[nh+1,tb] - soc_liion[1,tb+1] <= M * δ_inv_liion[tb]
    [tb in 1:ntb-1], soc_liion[1,tb+1] - liion.α_soc_max * r_liion[tb] <= M * (1 - δ_inv_liion[tb])
    [tb in 1:ntb-1], liion.α_soc_max * r_liion[tb] - soc_liion[1,tb+1] <= M * (1 - δ_inv_liion[tb])
    # SOH_liion
    [tb in 1:ntb-1], soh_liion[1,tb+1] - soh_liion[nh+1,tb] <= M * δ_inv_liion[tb]
    [tb in 1:ntb-1], soh_liion[nh+1,tb] - soh_liion[1,tb+1] <= M * δ_inv_liion[tb]
    [tb in 1:ntb-1], soh_liion[1,tb+1] - r_liion[tb] <= M * (1 - δ_inv_liion[tb])
    [tb in 1:ntb-1], r_liion[tb] - soh_liion[1,tb+1] <= M * (1 - δ_inv_liion[tb])

    # PV
    [tb in 1:ntb], pMax_pv[tb+1] - pMax_pv[tb] <= M * δ_inv_pv[tb]
    [tb in 1:ntb], pMax_pv[tb] - pMax_pv[tb+1] <= M * δ_inv_pv[tb]
    [tb in 1:ntb], pMax_pv[tb+1] - r_pv[tb] <= M * (1 - δ_inv_pv[tb])
    [tb in 1:ntb], r_pv[tb] - pMax_pv[tb+1] <= M * (1 - δ_inv_pv[tb])
    end)

    # Initial and final conditions
    @constraints(m, begin
    E_liion[1] == liion.Erated[1]
    pMax_pv[1] == pv.powerMax[1]
    soh_liion[1,1] == liion.soh[1,1] * liion.Erated[1]
    soc_liion[1,1] == liion.soc[1,1] * liion.Erated[1]
    [tb in 1:ntb], soc_liion[1,tb] == soc_liion[nh+1,tb]
    end)

    # Grid constraints
    @constraints(m, begin
    [h in 1:nh, tb in 2:ntb], p_g_in[h,tb] <= (1. - grid.τ_power) * maximum(ω_tb.clusters.ld_E[:,tb])
    [tb in 2:ntb], sum(p_g_in[h,tb] for h in 1:nh) <= (1. - grid.τ_energy) * sum(ω_tb.clusters.ld_E[h,tb] for h in 1:nh)
    end)

    # Objective
    @objective(m, Min, sum( γ[ω_tb.σ[tb]] * ( ω_tb.clusters.C_liion[tb] * r_liion[tb] + ω_tb.clusters.C_pv[tb] * r_pv[tb] +
    ω_tb.nby[tb] * sum( (p_g_in[h, tb] * ω_tb.clusters.C_grid_in[h,tb] + p_g_out[h, tb] * ω_tb.clusters.C_grid_out[h,tb]) * parameters.Δh for h=1:nh) ) for tb=1:ntb))

    return m
end
# With clustering
function multistage_milp_model(ld::Load, pv::Source, liion::Liion,
     controller::AnticipativeController, designer::AnticipativeMultiStageDesigner,
      grid::Grid, ω_tb::ClusteredScenarios, ω_td::ClusteredScenarios, parameters::NamedTuple)
    # Parameters
    M = 1000 # big-M value
    γ = 1. ./ (1. + parameters.τ) .^ (1:parameters.Y) # discount factor

    # Sets
    nh = size(ω_td.clusters.ld_E,1) # number of hours by td (=24)
    ntd = size(ω_td.clusters.ld_E,2) # number of td
    nd = size(ω_td.σ,1) # number of days over the year (length of the sequence)
    ntb = size(ω_td.clusters.ld_E,3) # number of tb

    # Model definition
    m = Model(CPLEX.Optimizer)

    # Variables
    # Operation control variables
    @variables(m, begin
    # Liion
    p_liion_ch[1:nh, 1:ntd, 1:ntb] <= 0.
    p_liion_dch[1:nh, 1:ntd, 1:ntb] >= 0.
    # Recourse
    p_g_out[1:nh, 1:ntd, 1:ntb] <= 0.
    p_g_in[1:nh, 1:ntd, 1:ntb] >= 0.
    end)
    # Investment control variables
    @variables(m, begin
    r_liion[1:ntb] >= 0.
    δ_inv_liion[1:ntb], Bin
    r_pv[1:ntb] >= 0.
    δ_inv_pv[1:ntb], Bin
    end)
    # Operation state variables
    @variables(m, begin
    soc_liion[1:nh+1, 1:nd, 1:ntb]
    soh_liion[1:nh+1, 1:nd, 1:ntb] >= 0.
    end)
    # Investment state variables
    @variables(m, begin
    E_liion[1:ntb+1] >= 0.
    pMax_pv[1:ntb+1] >= 0.
    end)
    # Operation constraints bounds
    @constraints(m, begin
    # Controls
    # Liion
    [h in 1:nh, td in 1:ntd, tb in 1:ntb], p_liion_dch[h,td,tb] <= liion.α_p_dch * E_liion[tb]
    [h in 1:nh, td in 1:ntd, tb in 1:ntb], p_liion_ch[h,td,tb] >= -liion.α_p_ch * E_liion[tb]
    # State
    # Liion
    [h in 1:nh+1, d in 1:nd, tb in 1:ntb], soc_liion[h,d,tb] <= liion.α_soc_max * E_liion[tb]
    [h in 1:nh+1, d in 1:nd, tb in 1:ntb], soc_liion[h,d,tb] >= liion.α_soc_min * E_liion[tb]
    [h in 1:nh+1, d in 1:nd, tb in 1:ntb], soh_liion[h,d,tb] <= E_liion[tb]
    end)
    # Investment constraints bounds
    @constraints(m, begin
    # Controls
    [tb in 1:ntb], r_liion[tb] <= 1000. * δ_inv_liion[tb]
    [tb in 1:ntb], r_pv[tb] <= 1000. * δ_inv_pv[tb]
    end)
    # Operation constraints dynamics and recourse
    @constraints(m, begin
    # Dynamics
    # Liion
    [h in 1:nh, d in 1:nd, tb in 1:ntb], soc_liion[h+1,d,tb] == soc_liion[h,d,tb] * (1 - liion.η_self * parameters.Δh) - (p_liion_ch[h,ω_td.σ[d,tb],tb] * liion.η_ch + p_liion_dch[h,ω_td.σ[d,tb],tb] / liion.η_dch) * parameters.Δh
    [d in 1:nd-1, tb in 1:ntb], soc_liion[1,d+1,tb] == soc_liion[nh+1,d,tb] * (1 - liion.η_self * parameters.Δh) - (p_liion_ch[1,ω_td.σ[d+1,tb],tb] * liion.η_ch + p_liion_dch[1,ω_td.σ[d+1,tb],tb] / liion.η_dch) * parameters.Δh
    [h in 1:nh, d in 1:nd, tb in 1:ntb], soh_liion[h+1,d,tb] == soh_liion[h,d,tb] - ω_tb.nby[tb] * (p_liion_dch[h,ω_td.σ[d,tb],tb] - p_liion_ch[h,ω_td.σ[d,tb],tb]) * parameters.Δh / (2. * liion.dod * liion.nCycle)
    [d in 1:nd-1, tb in 1:ntb], soh_liion[1,d+1,tb] == soh_liion[nh+1,d,tb] - ω_tb.nby[tb] * (p_liion_dch[1,ω_td.σ[d+1,tb],tb] - p_liion_ch[1,ω_td.σ[d+1,tb],tb]) * parameters.Δh / (2. * liion.dod * liion.nCycle)
    # Recourse
    [h in 1:nh, td in 1:ntd, tb in 1:ntb], ω_td.clusters.ld_E[h,td,tb] - pMax_pv[tb] * ω_td.clusters.pv_E[h,td,tb] <= p_g_out[h,td,tb] + p_g_in[h,td,tb] + p_liion_ch[h,td,tb] + p_liion_dch[h,td,tb]
    end)
    # Investment constraints dynamics
    @constraints(m, begin
    # Dynamics
    # E_liion
    # if δ_inv_liion = 0 (r_liion = 0) then E_liion[tb+1] =  E_liion[tb]
    [tb in 1:ntb], E_liion[tb+1] - E_liion[tb] <= M * δ_inv_liion[tb]
    [tb in 1:ntb], E_liion[tb] - E_liion[tb+1] <= M * δ_inv_liion[tb]
    # else E_liion[tb+1] = r_liion[tb] end
    [tb in 1:ntb], E_liion[tb+1] - r_liion[tb] <= M * (1 - δ_inv_liion[tb])
    [tb in 1:ntb], r_liion[tb] - E_liion[tb+1] <= M * (1 - δ_inv_liion[tb])
    # SOC_liion
    [tb in 1:ntb-1], soc_liion[1,1,tb+1] - soc_liion[nh+1,nd,tb] <= M * δ_inv_liion[tb]
    [tb in 1:ntb-1], soc_liion[nh+1,nd,tb] - soc_liion[1,1,tb+1] <= M * δ_inv_liion[tb]
    [tb in 1:ntb-1], soc_liion[1,1,tb+1] - liion.α_soc_max * r_liion[tb] <= M * (1 - δ_inv_liion[tb])
    [tb in 1:ntb-1], liion.α_soc_max * r_liion[tb] - soc_liion[1,1,tb+1] <= M * (1 - δ_inv_liion[tb])
    # SOH_liion
    [tb in 1:ntb-1], soh_liion[1,1,tb+1] - soh_liion[nh+1,nd,tb] <= M * δ_inv_liion[tb]
    [tb in 1:ntb-1], soh_liion[nh+1,nd,tb] - soh_liion[1,1,tb+1] <= M * δ_inv_liion[tb]
    [tb in 1:ntb-1], soh_liion[1,1,tb+1] - r_liion[tb] <= M * (1 - δ_inv_liion[tb])
    [tb in 1:ntb-1], r_liion[tb] - soh_liion[1,1,tb+1] <= M * (1 - δ_inv_liion[tb])
    # PV
    [tb in 1:ntb], pMax_pv[tb+1] - pMax_pv[tb] <= M * δ_inv_pv[tb]
    [tb in 1:ntb], pMax_pv[tb] - pMax_pv[tb+1] <= M * δ_inv_pv[tb]
    [tb in 1:ntb], pMax_pv[tb+1] - r_pv[tb] <= M * (1 - δ_inv_pv[tb])
    [tb in 1:ntb], r_pv[tb] - pMax_pv[tb+1] <= M * (1 - δ_inv_pv[tb])
    end)
    # Initial and final conditions
    @constraints(m, begin
    E_liion[1] == liion.Erated[1]
    pMax_pv[1] == pv.powerMax[1]
    soh_liion[1,1,1] == liion.soh[1,1] * liion.Erated[1]
    soc_liion[1,1,1] == liion.soc[1,1] * liion.Erated[1]
    [tb in 1:ntb], soc_liion[1,1,tb] == soc_liion[nh+1,nd,tb]
    end)
    # Grid constraints
    @constraints(m, begin
    [h in 1:nh, td in 1:ntd, tb in 2:ntb], p_g_in[h,td,tb] <= (1. - grid.τ_power) * maximum(ω_tb.clusters.ld_E[:,tb])
    [tb in 2:ntb], sum(ω_td.nby[td,tb] * sum(p_g_in[h,td,tb] for h in 1:nh) for td in 1:ntd) <= (1. - grid.τ_energy) * sum(ω_tb.clusters.ld_E[:,tb])
    end)

    # Objective
    @objective(m, Min, sum( γ[ω_tb.σ[tb]] * ( ω_tb.clusters.C_liion[tb] * r_liion[tb] + ω_tb.clusters.C_pv[tb] * r_pv[tb] +
    ω_tb.nby[tb] * sum(ω_td.nby[td,tb] * sum((p_g_in[h,td,tb] * ω_tb.clusters.C_grid_in[h,tb] + p_g_out[h,td,tb] * ω_tb.clusters.C_grid_out[h,tb]) * parameters.Δh for h=1:nh) for td in 1:ntd) ) for tb=1:ntb))

    return m
end

#                                   Offline functions
#_______________________________________________________________________________
function offline_optimization(ld::Load, pv::Source, liion::Liion,
    controller::AnticipativeController, designer::AnticipativeMultiStageDesigner,
     grid::Grid, ω_optim::Scenarios, ntb::Int64, parameters::NamedTuple)

     # Selection of one-year scenarios from the optimization dataset
     ω_milp = Scenarios(ω_optim.timestamp,ω_optim.ld_E[:,:,1],ω_optim.ld_H[:,:,1],ω_optim.pv_E[:,:,1],
     ω_optim.C_pv[:,1],ω_optim.C_liion[:,1],ω_optim.C_tes[:,1],ω_optim.C_tank[:,1],
     ω_optim.C_elyz[:,1],ω_optim.C_fc[:,1],ω_optim.C_heater[:,1],
     ω_optim.C_grid_in[:,:,1],ω_optim.C_grid_out[:,:,1])

     # Time block clustering
     ω_tb = clustering_time_block(ω_milp, ntb)

     # Initialize model
     designer.model = controller.model = multistage_milp_model(
     ld, pv, liion, controller, designer, grid, ω_tb, parameters)

     # Compute both investment and operation decisions
     optimize!(designer.model)

     # Formatting variables to simulation
     # Operation decisions
     controller.u = (
     u_liion =  reshape_tb_to_simu_operation(value.(controller.model[:p_liion_ch] + controller.model[:p_liion_dch]), ω_tb.nby),
     )
     # Investment decisions
     designer.u =(
     u_liion = reshape_tb_to_simu_investment(value.(designer.model[:r_liion]), ω_tb.σ),
     u_pv = reshape_tb_to_simu_investment(value.(designer.model[:r_pv]), ω_tb.σ),
     )
end
# With clustering
function offline_optimization(ld::Load, pv::Source, liion::Liion,
    controller::AnticipativeController, designer::AnticipativeMultiStageDesigner,
     grid::Grid, ω_optim::Scenarios, ntb::Int64, ntd::Int64, parameters::NamedTuple)

     # Selection of one-year scenarios from the optimization dataset
     ω_milp = Scenarios(ω_optim.timestamp,ω_optim.ld_E[:,:,1],nothing,ω_optim.pv_E[:,:,1],
     ω_optim.C_pv[:,1],ω_optim.C_liion[:,1],nothing,nothing,nothing,nothing,nothing,
     ω_optim.C_grid_in[:,:,1],ω_optim.C_grid_out[:,:,1])

     # Time block clustering
     ω_tb = clustering_time_block(ω_milp, ntb)

     # Typical days clustering
     ω_td = clustering_typical_day(ω_tb, ntd)

     # Initialize investment model
     designer.model = multistage_milp_model(ld, pv, liion, controller, designer,
      grid, ω_tb, ω_td, parameters)

     # Compute investment decisions with typical days
     optimize!(designer.model)

     # Initialize operation model with previous sizes
     controller.model = multistage_milp_model(ld, pv, liion, controller, designer,
      grid, ω_tb, parameters)

     # TODO: pas de convergence, peut etre enlever contrainte réseau ?
     # Fix sizes from investment optimization
     fix.(controller.model[:r_liion], value.(designer.model[:r_liion]), force = true)
     fix.(controller.model[:δ_inv_liion], value.(designer.model[:δ_inv_liion]), force = true)
     fix.(controller.model[:r_pv], value.(designer.model[:r_pv]), force = true)
     fix.(controller.model[:δ_inv_pv], value.(designer.model[:δ_inv_pv]), force = true)

     # Compute operation decisions with initial scenarios
     optimize!(controller.model)

     # Formatting variables to simulation
     # Operation decisions
     controller.u = (
     u_liion =  reshape_tb_to_simu_operation(value.(controller.model[:p_liion_ch] + controller.model[:p_liion_dch]), ω_tb.nby),
     )
     # Investment decisions
     designer.u =(
     u_liion = reshape_tb_to_simu_investment(value.(designer.model[:r_liion]), ω_tb.σ),
     u_pv = reshape_tb_to_simu_investment(value.(designer.model[:r_pv]), ω_tb.σ),
     )
end

#                                   Online functions
#_______________________________________________________________________________
function compute_investment_decisions(y::Int64, s::Int64, ld::Load, pv::Source,
    liion::Liion, grid::Grid, controller::AnticipativeController,
    designer::AnticipativeMultiStageDesigner, ω_optim::Scenarios, parameters::NamedTuple)
    return nothing
end
