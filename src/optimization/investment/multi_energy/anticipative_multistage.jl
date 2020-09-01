#=
    Anticipative multi-stage designer with anticipative controller
=#

# TODO : Mettre à jour ou supprimer !

#                                   Model definition
#_______________________________________________________________________________
function multistage_milp_model(ld::Load, pv::Source, liion::Liion, h2tank::H2Tank,
   elyz::Electrolyzer, fc::FuelCell, tes::ThermalSto, heater::Heater,
    controller::AnticipativeController, designer::AnticipativeMultiStageDesigner,
     grid:: Grid, γ)

    # Parameters
    ϵ = 0.001
    M1 = 1000
    M2 = 20000
    M3 = 50

    # Sets
    # Number of time blocs
    nyy = length(1:designer.Δy:designer.horizon) + 1
    # Number of operation variables
    nh = length(1:controller.Δh:controller.horizon)

    # Sequence of decisions
    σ_u = vcat(collect(1:designer.Δy:designer.horizon),designer.horizon)

    # Model definition
    m = Model(with_optimizer(CPLEX.Optimizer))
    MOI.set(m, MOI.RawParameter("CPXPARAM_MIP_Cuts_Gomory"), 2)
    MOI.set(m, MOI.RawParameter("CPXPARAM_MIP_Strategy_Probe"), 3)
    MOI.set(m, MOI.RawParameter("CPXPARAM_MIP_Tolerances_MIPGap"), 0.05)

    # Variables
    # Operation control variables
    @variables(m, begin
    # Liion
     p_liion_ch[1:nh, 1:nyy] <= 0.
     p_liion_dch[1:nh, 1:nyy] >= 0.
     # H2 hub
     p_elyz_E[1:nh, 1:nyy] <= 0.
     δ_elyz[1:nh, 1:nyy], (Bin)
     p_fc_E[1:nh, 1:nyy] >= 0.
     δ_fc[1:nh, 1:nyy], (Bin)
     # TES
     p_tes_ch[1:nh, 1:nyy] <= 0.
     p_tes_dch[1:nh, 1:nyy] >= 0.
     # Heater
     p_heater_E[1:nh, 1:nyy] <= 0.
     # Recourse
     p_g_out[1:nh, 1:nyy] <= 0.
     p_g_in[1:nh, 1:nyy] >= 0.
    end)
    # Investment control variables
    @variables(m, begin
    r_pv[1:nyy] >= 0.
    δ_inv_pv[1:nyy], (Bin)
    r_liion[1:nyy] >= 0.
    δ_inv_liion[1:nyy], (Bin)
    r_tank[1:nyy] >= 0.
    δ_inv_tank[1:nyy], (Bin)
    r_elyz[1:nyy] >= 0.
    δ_inv_elyz[1:nyy], (Bin)
    r_fc[1:nyy] >= 0.
    δ_inv_fc[1:nyy], (Bin)
    r_tes[1:nyy] >= 0.
    δ_inv_tes[1:nyy], (Bin)
    end)
    # Operation state variables
    @variables(m, begin
    soc_liion[1:nh+1, 1:nyy]
    soc_h2[1:nh+1, 1:nyy]
    soc_tes[1:nh+1, 1:nyy]
    soh_liion[1:nh+1, 1:nyy] >= 0.
    soh_elyz[1:nh+1, 1:nyy] >= 0.
    soh_fc[1:nh+1, 1:nyy] >= 0.
    end)
    # Investment state variables
    @variables(m, begin
    pMax_pv[1:nyy+1] >= 0.
    E_liion[1:nyy+1] >= 0.
    E_h2[1:nyy+1] >= 0.
    pMax_elyz[1:nyy+1] >= 0.
    pMax_fc[1:nyy+1] >= 0.
    E_tes[1:nyy+1] >= 0.
    end)

    # Operation constraints bounds
    @constraints(m, begin
    # Controls
    # Storage
    # Liion
    [h in 1:nh, y in 1:nyy], p_liion_dch[h,y] <= liion.α_p_dch * E_liion[y]
    [h in 1:nh, y in 1:nyy], p_liion_ch[h,y] >= -liion.α_p_ch * E_liion[y]
    # TES
    [h in 1:nh, y in 1:nyy], p_tes_dch[h,y] <= tes.α_p_dch * E_tes[y]
    [h in 1:nh, y in 1:nyy], p_tes_ch[h,y] >= -tes.α_p_ch * E_tes[y]

    # Converters
    [h in 1:nh, y in 1:nyy], p_fc_E[h,y] <= pMax_fc[y]
    [h in 1:nh, y in 1:nyy], p_fc_E[h,y] <= ϵ + M3 * δ_fc[h,y]
    [h in 1:nh, y in 1:nyy], p_elyz_E[h,y] >= -pMax_elyz[y]
    [h in 1:nh, y in 1:nyy], p_elyz_E[h,y] >= -ϵ - M3 * δ_elyz[h,y]
    [h in 1:nh, y in 1:nyy], p_heater_E[h,y] >= -heater.powerMax[1]


    # State
    # Liion
    [h in 1:nh+1, y in 1:nyy], soc_liion[h,y] <= liion.α_soc_max * E_liion[y]
    [h in 1:nh+1, y in 1:nyy], soc_liion[h,y] >= liion.α_soc_min * E_liion[y]
    [h in 1:nh+1, y in 1:nyy], soh_liion[h,y] <= E_liion[y]
    # H2 hub
    [h in 1:nh+1, y in 1:nyy], soc_h2[h,y] <= h2tank.α_soc_max * E_h2[y]
    [h in 1:nh+1, y in 1:nyy], soc_h2[h,y] >= h2tank.α_soc_min * E_h2[y]
    [h in 1:nh+1, y in 1:nyy], soh_elyz[h,y] <= 1.
    [h in 1:nh+1, y in 1:nyy], soh_fc[h,y] <= 1.
    # TES
    [h in 1:nh+1, y in 1:nyy], soc_tes[h,y] <= tes.α_soc_max * E_tes[y]
    [h in 1:nh+1, y in 1:nyy], soc_tes[h,y] >= tes.α_soc_min * E_tes[y]
    end)
    # Investment constraints bounds
    @constraints(m, begin
    # Controls
    [y in 1:nyy], r_pv[y] <= 1000 * δ_inv_pv[y]
    [y in 1:nyy], r_liion[y] <= 1000 * δ_inv_liion[y]
    [y in 1:nyy], r_tank[y] <= 20000 * δ_inv_tank[y]
    [y in 1:nyy], r_elyz[y] <= 50 * δ_inv_elyz[y]
    [y in 1:nyy], r_fc[y] <= 50 * δ_inv_fc[y]
    [y in 1:nyy], r_tes[y] <= 1000 * δ_inv_tes[y]
    end)
    # Operation constraints dynamics and recourse
    @constraints(m, begin
    # Dynamics
    # Liion
    [h in 1:nh, y in 1:nyy], soc_liion[h+1,y] == soc_liion[h,y] * (1 - liion.η_self * controller.Δh) - (p_liion_ch[h,y] * liion.η_ch + p_liion_dch[h,y] / liion.η_dch) * controller.Δh
    [h in 1:nh, y in 1:nyy], soh_liion[h+1,y] == soh_liion[h,y] - designer.Δy  * (p_liion_dch[h,y] - p_liion_ch[h,y]) * controller.Δh / (2. * liion.dod * liion.nCycle)
    # H2 hub
    [h in 1:nh, y in 1:nyy], soc_h2[h+1,y] == soc_h2[h,y] * (1 - h2tank.η_self * controller.Δh)  - (p_elyz_E[h,y] * elyz.η_E_H2 + p_fc_E[h,y] / fc.η_H2_E) * controller.Δh
    [h in 1:nh, y in 1:nyy], soh_elyz[h+1,y] == soh_elyz[h,y] - designer.Δy  * δ_elyz[h,y] * controller.Δh / elyz.nHoursMax
    [h in 1:nh, y in 1:nyy], soh_fc[h+1,y] == soh_fc[h,y] - designer.Δy  * δ_fc[h,y] * controller.Δh / fc.nHoursMax
    # TES
    [h in 1:nh, y in 1:nyy], soc_tes[h+1,y] == soc_tes[h,y] * (1 - tes.η_self * controller.Δh)  - (p_tes_ch[h,y] * tes.η_ch + p_tes_dch[h,y] / tes.η_dch) * controller.Δh

    # Recourse
    [h in 1:nh, y in 1:nyy], ld.power_E[h,σ_u[y]] - pMax_pv[y] * pv.power_E[h,σ_u[y]] <= p_g_out[h,y] + p_g_in[h,y] + p_liion_ch[h,y] + p_liion_dch[h,y] + p_elyz_E[h,y] + p_fc_E[h,y] + p_heater_E[h,y]
    [h in 1:nh, y in 1:nyy], ld.power_H[h,σ_u[y]] <= - elyz.η_E_H * p_elyz_E[h,y] +  fc.η_H2_H / fc.η_H2_E * p_fc_E[h,y] -  heater.η_E_H * p_heater_E[h,y] + p_tes_ch[h,y] + p_tes_dch[h,y]
    end)
    # Investment constraints dynamics
    @constraints(m, begin
    # Dynamics
    # PV
    [y in 1:nyy], pMax_pv[y+1] - pMax_pv[y] <= M1 * δ_inv_pv[y]
    [y in 1:nyy], pMax_pv[y] - pMax_pv[y+1] <= M1 * δ_inv_pv[y]
    [y in 1:nyy], pMax_pv[y+1] - r_pv[y] <= M1 * (1 - δ_inv_pv[y])
    [y in 1:nyy], r_pv[y] - pMax_pv[y+1] <= M1 * (1 - δ_inv_pv[y])

    # Liion
    # E_liion
    # if δ_inv_liion = 0 (r_liion = 0) then E_liion[y+1] =  E_liion[y]
    [y in 1:nyy], E_liion[y+1] - E_liion[y] <= M1 * δ_inv_liion[y]
    [y in 1:nyy], E_liion[y] - E_liion[y+1] <= M1 * δ_inv_liion[y]
    # else E_liion[y+1] = r_liion[y] end
    [y in 1:nyy], E_liion[y+1] - r_liion[y] <= M1 * (1 - δ_inv_liion[y])
    [y in 1:nyy], r_liion[y] - E_liion[y+1] <= M1 * (1 - δ_inv_liion[y])
    # SOC_liion
    [y in 1:nyy-1], soc_liion[1,y+1] - soc_liion[nh+1,y] <= M1 * δ_inv_liion[y]
    [y in 1:nyy-1], soc_liion[nh+1,y] - soc_liion[1,y+1] <= M1 * δ_inv_liion[y]
    [y in 1:nyy-1], soc_liion[1,y+1] - liion.α_soc_max * r_liion[y] <= M1 * (1 - δ_inv_liion[y])
    [y in 1:nyy-1], liion.α_soc_max * r_liion[y] - soc_liion[1,y+1] <= M1 * (1 - δ_inv_liion[y])
    # SOH_liion
    [y in 1:nyy-1], soh_liion[1,y+1] - soh_liion[nh+1,y] <= M1 * δ_inv_liion[y]
    [y in 1:nyy-1], soh_liion[nh+1,y] - soh_liion[1,y+1] <= M1 * δ_inv_liion[y]
    [y in 1:nyy-1], soh_liion[1,y+1] - r_liion[y] <= M1 * (1 - δ_inv_liion[y])
    [y in 1:nyy-1], r_liion[y] - soh_liion[1,y+1] <= M1 * (1 - δ_inv_liion[y])

    # H2 hub
    # E_h2
    [y in 1:nyy], E_h2[y+1] - E_h2[y] <= M2 * δ_inv_tank[y]
    [y in 1:nyy], E_h2[y] - E_h2[y+1] <= M2 * δ_inv_tank[y]
    [y in 1:nyy], E_h2[y+1] - r_tank[y] <= M2 * (1 - δ_inv_tank[y])
    [y in 1:nyy], r_tank[y] - E_h2[y+1] <= M2 * (1 - δ_inv_tank[y])
    # SOC_h2
    [y in 1:nyy-1], soc_h2[1,y+1] - soc_h2[nh+1,y] <= M2 * δ_inv_tank[y]
    [y in 1:nyy-1], soc_h2[nh+1,y] - soc_h2[1,y+1] <= M2 * δ_inv_tank[y]
    [y in 1:nyy-1], soc_h2[1,y+1] - h2tank.α_soc_max * r_tank[y] <= M2 * (1 - δ_inv_tank[y])
    [y in 1:nyy-1], h2tank.α_soc_max * r_tank[y] - soc_h2[1,y+1] <= M2 * (1 - δ_inv_tank[y])
    # pMax_elyz
    [y in 1:nyy], pMax_elyz[y+1] - pMax_elyz[y] <= M3 * δ_inv_elyz[y]
    [y in 1:nyy], pMax_elyz[y] - pMax_elyz[y+1] <= M3 * δ_inv_elyz[y]
    [y in 1:nyy], pMax_elyz[y+1] - r_elyz[y] <= M3 * (1 - δ_inv_elyz[y])
    [y in 1:nyy], r_elyz[y] - pMax_elyz[y+1] <= M3 * (1 - δ_inv_elyz[y])
    # SOH_elyz
    [y in 1:nyy-1], soh_elyz[1,y+1] - soh_elyz[nh+1,y] <= M3 * δ_inv_elyz[y]
    [y in 1:nyy-1], soh_elyz[nh+1,y] - soh_elyz[1,y+1] <= M3 * δ_inv_elyz[y]
    [y in 1:nyy-1], soh_elyz[1,y+1] - 1. <= M3 * (1 - δ_inv_elyz[y])
    [y in 1:nyy-1], 1. - soh_elyz[1,y+1] <= M3 * (1 - δ_inv_elyz[y])
    # pMax_fc
    [y in 1:nyy], pMax_fc[y+1] - pMax_fc[y] <= M3 * δ_inv_fc[y]
    [y in 1:nyy], pMax_fc[y] - pMax_fc[y+1] <= M3 * δ_inv_fc[y]
    [y in 1:nyy], pMax_fc[y+1] - r_fc[y] <= M3 * (1 - δ_inv_fc[y])
    [y in 1:nyy], r_fc[y] - pMax_fc[y+1] <= M3 * (1 - δ_inv_fc[y])
    # SOH_fc
    [y in 1:nyy-1], soh_fc[1,y+1] - soh_fc[nh+1,y] <= M3 * δ_inv_fc[y]
    [y in 1:nyy-1], soh_fc[nh+1,y] - soh_fc[1,y+1] <= M3 * δ_inv_fc[y]
    [y in 1:nyy-1], soh_fc[1,y+1] - 1. <= M3 * (1 - δ_inv_fc[y])
    [y in 1:nyy-1], 1. - soh_fc[1,y+1] <= M3 * (1 - δ_inv_fc[y])

    # TES
    # E_tes
    [y in 1:nyy], E_tes[y+1] - E_tes[y] <= M1 * δ_inv_tes[y]
    [y in 1:nyy], E_tes[y] - E_tes[y+1] <= M1 * δ_inv_tes[y]
    [y in 1:nyy], E_tes[y+1] - r_tes[y] <= M1 * (1 - δ_inv_tes[y])
    [y in 1:nyy], r_tes[y] - E_tes[y+1] <= M1 * (1 - δ_inv_tes[y])
    # SOC_tes
    [y in 1:nyy-1], soc_tes[1,y+1] - soc_tes[nh+1,y] <= M1 * δ_inv_tes[y]
    [y in 1:nyy-1], soc_tes[nh+1,y] - soc_tes[1,y+1] <= M1 * δ_inv_tes[y]
    [y in 1:nyy-1], soc_tes[1,y+1] - tes.α_soc_max * r_tes[y] <= M1 * (1 - δ_inv_tes[y])
    [y in 1:nyy-1], tes.α_soc_max * r_tes[y] - soc_tes[1,y+1] <= M1 * (1 - δ_inv_tes[y])
    end)

    # Initial and final conditions
    @constraints(m, begin
    # Liion
    E_liion[1] == liion.Erated[1]
    soh_liion[1,1] == liion.soh[1,1] * liion.Erated[1]
    soc_liion[1,1] == liion.soc[1,1] * liion.Erated[1]
    # [y in 1:nyy], soc_liion[1,y] == soc_liion[nh+1,y]
    # H2
    E_h2[1] == h2tank.Erated[1]
    soc_h2[1,1] == h2tank.soc[1,1] * h2tank.Erated[1]
    # [y in 1:nyy], soc_h2[1,y] == soc_h2[nh+1,y]
    # TES
    E_tes[1] == tes.Erated[1]
    soc_tes[1,1] == tes.soc[1,1] * tes.Erated[1]
    # [y in 1:nyy], soc_tes[1,y] == soc_tes[nh+1,y]
    # PV
    pMax_pv[1] == pv.powerMax[1]
    # Elyz
    pMax_elyz[1] == elyz.powerMax[1]
    soh_elyz[1,1] == elyz.soh[1,1]
    # FC
    pMax_fc[1] == fc.powerMax[1]
    soh_fc[1,1] == fc.soh[1,1]

    end)

    # Grid constraints
    @constraints(m, begin
    [h in 1:nh, y in 2:nyy], p_g_in[h, y] <= (1. - grid.τ_power) * maximum(ld.power_E[:,σ_u[y]] .+ ld.power_H[:,σ_u[y]] / heater.η_E_H)
    [y in 2:nyy], sum(p_g_in[h,y] for h in 1:nh) <= (1. - grid.τ_energy) * sum(ld.power_E[h,σ_u[y]] for h in 1:nh)
    end)

    # Objective
    @objective(m, Min, sum( γ[σ_u[y]] * (pv.C_pv[σ_u[y]] * r_pv[y] + liion.C_liion[σ_u[y]] * r_liion[y] + h2tank.C_tank[σ_u[y]] * r_tank[y] +
    elyz.C_elyz[σ_u[y]] * r_elyz[y] + fc.C_fc[σ_u[y]] * r_fc[y] + tes.C_tes[σ_u[y]] * r_tes[y] +
    designer.Δy * sum( (p_g_in[h, y] * grid.C_grid_in[h, σ_u[y]] + p_g_out[h, y] * grid.C_grid_out[h, σ_u[y]]) * controller.Δh  for h=1:nh) ) for y=1:nyy))

    # Solve
    optimize!(m)

    # Outputs
    # Operation controls
    controller.u = (
    u_liion =  reshape(repeat( value.(p_liion_ch + p_liion_dch), ceil(Int,nh/nh) * designer.Δy,1),nh,:)[:,designer.Δy:end],
    u_elyz = reshape(repeat( value.(p_elyz_E), ceil(Int,nh/nh) * designer.Δy,1),nh,:)[:,designer.Δy:end],
    u_fc= reshape(repeat( value.(p_fc_E), ceil(Int,nh/nh) * designer.Δy,1),nh,:)[:,designer.Δy:end],
    u_tes = reshape(repeat( value.(p_tes_ch + p_tes_dch), ceil(Int,nh/nh) * designer.Δy,1),nh,:)[:,designer.Δy:end],
    u_heater = reshape(repeat( value.(p_heater_E), ceil(Int,nh/nh) * designer.Δy,1),nh,:)[:,designer.Δy:end],
    )
    # Investment controls
    designer.u = (
    u_pv = reshape(vcat(zeros(designer.Δy-1,nyy),reshape(value.(r_pv),1,:)),:,1)[designer.Δy:end],
    u_liion = reshape(vcat(zeros(designer.Δy-1,nyy),reshape(value.(r_liion),1,:)),:,1)[designer.Δy:end],
    u_tank = reshape(vcat(zeros(designer.Δy-1,nyy),reshape(value.(r_tank),1,:)),:,1)[designer.Δy:end],
    u_elyz = reshape(vcat(zeros(designer.Δy-1,nyy),reshape(value.(r_elyz),1,:)),:,1)[designer.Δy:end],
    u_fc = reshape(vcat(zeros(designer.Δy-1,nyy),reshape(value.(r_fc),1,:)),:,1)[designer.Δy:end],
    u_tes = reshape(vcat(zeros(designer.Δy-1,nyy),reshape(value.(r_tes),1,:)),:,1)[designer.Δy:end],
    )
    # Operation states
    controller.x = (
    soc_liion = value.(soc_liion),
    soc_h2 = value.(soc_h2),
    soh_elyz = value.(soh_elyz),
    soh_fc = value.(soh_fc),
    soc_tes = value.(soc_tes),
    soh_liion = value.(soh_liion),
    )
    # Investment states
    designer.x = (
    pMax_pv = value.(pMax_pv),
    E_liion = value.(E_liion),
    E_h2 = value.(E_h2),
    pMax_elyz = value.(pMax_elyz),
    pMax_fc = value.(pMax_fc),
    E_tes = value.(E_tes),
    )

    # Log
    controller.log = designer.log = termination_status(m)

    # Model
    controller.model = designer.model = m
end
# With clustering
function multistage_milp_model(ld::Load, pv::Source, liion::Liion, h2tank::H2Tank,
   elyz::Electrolyzer, fc::FuelCell, tes::ThermalSto, heater::Heater,
    controller::AnticipativeController, designer::AnticipativeMultiStageDesigner,
     grid:: Grid, γ, nCluster)
    # Parameters
    ϵ = 0.001
    M1 = 1000
    M2 = 20000
    M3 = 50

    # CLustering
    ld_E_bytd, ld_H_bytd, pv_bytd, σ_td, n_bytd = clustering_data(ld, pv, nCluster)
    # Sets
    nh = size(ld_E_bytd,1) # number of hours by td (=24)
    ntd = size(ld_E_bytd,2) # number of td (=nCluster)
    nd = size(σ_td,1) # number of days over the year (length of the sequence)
    nyy = length(1:designer.Δy:designer.horizon) + 1 # number of decisions for the design
    # Sequence of design decisions
    σ_u = vcat(collect(1:designer.Δy:designer.horizon),designer.horizon)
    # Model definition
    m = Model(with_optimizer(CPLEX.Optimizer))
    # MOI.set(m, MOI.RawParameter("CPXPARAM_MIP_Cuts_Gomory"), 2)
    # MOI.set(m, MOI.RawParameter("CPXPARAM_MIP_Strategy_Probe"), 3)
    MOI.set(m, MOI.RawParameter("CPXPARAM_MIP_Tolerances_MIPGap"), 0.05)
    # Variables
    # Operation control variables
    @variables(m, begin
    # Liion
     p_liion_ch[1:nh, 1:ntd, 1:nyy] <= 0.
     p_liion_dch[1:nh, 1:ntd, 1:nyy] >= 0.
     # H2 hub
     p_elyz_E[1:nh, 1:ntd, 1:nyy] <= 0.
     δ_elyz[1:nh, 1:ntd, 1:nyy], (Bin)
     p_fc_E[1:nh, 1:ntd, 1:nyy] >= 0.
     δ_fc[1:nh, 1:ntd, 1:nyy], (Bin)
     # TES
     p_tes_ch[1:nh, 1:ntd, 1:nyy] <= 0.
     p_tes_dch[1:nh, 1:ntd, 1:nyy] >= 0.
     # Heater
     p_heater_E[1:nh, 1:ntd, 1:nyy] <= 0.
     # Recourse
     p_g_out[1:nh, 1:ntd, 1:nyy] <= 0.
     p_g_in[1:nh, 1:ntd, 1:nyy] >= 0.
    end)
    # Investment control variables
    @variables(m, begin
    r_pv[1:nyy] >= 0.
    δ_inv_pv[1:nyy], (Bin)
    r_liion[1:nyy] >= 0.
    δ_inv_liion[1:nyy], (Bin)
    r_tank[1:nyy] >= 0.
    δ_inv_tank[1:nyy], (Bin)
    r_elyz[1:nyy] >= 0.
    δ_inv_elyz[1:nyy], (Bin)
    r_fc[1:nyy] >= 0.
    δ_inv_fc[1:nyy], (Bin)
    r_tes[1:nyy] >= 0.
    δ_inv_tes[1:nyy], (Bin)
    end)
    # Operation state variables
    @variables(m, begin
    soc_liion[1:nh+1, 1:nd, 1:nyy]
    soc_h2[1:nh+1, 1:nd, 1:nyy]
    soc_tes[1:nh+1, 1:nd, 1:nyy]
    soh_liion[1:nh+1, 1:nd, 1:nyy] >= 0.
    soh_elyz[1:nh+1, 1:nd, 1:nyy] >= 0.
    soh_fc[1:nh+1, 1:nd, 1:nyy] >= 0.
    end)
    # Investment state variables
    @variables(m, begin
    pMax_pv[1:nyy+1] >= 0.
    E_liion[1:nyy+1] >= 0.
    E_h2[1:nyy+1] >= 0.
    pMax_elyz[1:nyy+1] >= 0.
    pMax_fc[1:nyy+1] >= 0.
    E_tes[1:nyy+1] >= 0.
    end)
    # Operation constraints bounds
    @constraints(m, begin
    # Controls
    # Storage
    # Liion
    [h in 1:nh, td in 1:ntd, y in 1:nyy], p_liion_dch[h,td,y] <= liion.α_p_dch * E_liion[y]
    [h in 1:nh, td in 1:ntd, y in 1:nyy], p_liion_ch[h,td,y] >= -liion.α_p_ch * E_liion[y]
    # TES
    [h in 1:nh, td in 1:ntd, y in 1:nyy], p_tes_dch[h,td,y] <= tes.α_p_dch * E_tes[y]
    [h in 1:nh, td in 1:ntd, y in 1:nyy], p_tes_ch[h,td,y] >= -tes.α_p_ch * E_tes[y]
    # Converters
    [h in 1:nh, td in 1:ntd, y in 1:nyy], p_fc_E[h,td,y] <= pMax_fc[y]
    [h in 1:nh, td in 1:ntd, y in 1:nyy], p_fc_E[h,td,y] <= ϵ + M3 * δ_fc[h,td,y]
    [h in 1:nh, td in 1:ntd, y in 1:nyy], p_elyz_E[h,td,y] >= -pMax_elyz[y]
    [h in 1:nh, td in 1:ntd, y in 1:nyy], p_elyz_E[h,td,y] >= -ϵ - M3 * δ_elyz[h,td,y]
    [h in 1:nh, td in 1:ntd, y in 1:nyy], p_heater_E[h,td,y] >= -heater.powerMax[1]
    # State
    # Liion
    [h in 1:nh+1, d in 1:nd, y in 1:nyy], soc_liion[h,d,y] <= liion.α_soc_max * E_liion[y]
    [h in 1:nh+1, d in 1:nd, y in 1:nyy], soc_liion[h,d,y] >=liion.α_soc_min * E_liion[y]
    [h in 1:nh+1, d in 1:nd, y in 1:nyy], soh_liion[h,d,y] <= E_liion[y]
    # H2 hub
    [h in 1:nh+1, d in 1:nd, y in 1:nyy], soc_h2[h,d,y] <= h2tank.α_soc_max * E_h2[y]
    [h in 1:nh+1, d in 1:nd, y in 1:nyy], soc_h2[h,d,y] >= h2tank.α_soc_min * E_h2[y]
    [h in 1:nh+1, d in 1:nd, y in 1:nyy], soh_elyz[h,d,y] <= 1.
    [h in 1:nh+1, d in 1:nd, y in 1:nyy], soh_fc[h,d,y] <= 1.
    # TES
    [h in 1:nh+1, d in 1:nd, y in 1:nyy], soc_tes[h,d,y] <= tes.α_soc_max * E_tes[y]
    [h in 1:nh+1, d in 1:nd, y in 1:nyy], soc_tes[h,d,y] >= tes.α_soc_min * E_tes[y]
    end)
    # Investment constraints bounds
    @constraints(m, begin
    # Controls
    [y in 1:nyy], r_pv[y] <= 1000 * δ_inv_pv[y]
    [y in 1:nyy], r_liion[y] <= 1000 * δ_inv_liion[y]
    [y in 1:nyy], r_tank[y] <= 20000 * δ_inv_tank[y]
    [y in 1:nyy], r_elyz[y] <= 50 * δ_inv_elyz[y]
    [y in 1:nyy], r_fc[y] <= 50 * δ_inv_fc[y]
    [y in 1:nyy], r_tes[y] <= 1000 * δ_inv_tes[y]
    end)
    # Operation constraints dynamics and recourse
    @constraints(m, begin
    # Dynamics
    # Liion dynamics by td
    [h in 1:nh, d in 1:nd, y in 1:nyy], soc_liion[h+1,d,y] == soc_liion[h,d,y] * (1 - liion.η_self * controller.Δh) - (p_liion_ch[h,σ_td[d,σ_u[y]],y] * liion.η_ch + p_liion_dch[h,σ_td[d,σ_u[y]],y] / liion.η_dch) * controller.Δh
    [h in 1:nh, d in 1:nd, y in 1:nyy], soh_liion[h+1,d,y] == soh_liion[h,d,y] - designer.Δy * (p_liion_dch[h,σ_td[d,σ_u[y]],y] - p_liion_ch[h,σ_td[d,σ_u[y]],y]) * controller.Δh / (2. * liion.dod * liion.nCycle)
    # Liion coupling between td
    [d in 1:nd-1, y in 1:nyy], soc_liion[1,d+1,y] == soc_liion[nh+1,d,y] * (1 - liion.η_self * controller.Δh) - (p_liion_ch[1,σ_td[d+1,σ_u[y]],y] * liion.η_ch + p_liion_dch[1,σ_td[d+1,σ_u[y]],y] / liion.η_dch) * controller.Δh
    [d in 1:nd-1, y in 1:nyy], soh_liion[1,d+1,y] == soh_liion[nh+1,d,y] - designer.Δy * (p_liion_dch[1,σ_td[d+1,σ_u[y]],y] - p_liion_ch[1,σ_td[d+1,σ_u[y]],y]) * controller.Δh / (2. * liion.dod * liion.nCycle)
    # H2 tank dynamic by td
    [h in 1:nh, d in 1:nd, y in 1:nyy], soc_h2[h+1,d,y] == soc_h2[h,d,y] * (1 - h2tank.η_self * controller.Δh)  - (p_elyz_E[h,σ_td[d,σ_u[y]],y] * elyz.η_E_H2 + p_fc_E[h,σ_td[d,σ_u[y]],y] / fc.η_H2_E) * controller.Δh
    # H2 tank coupling between td
    [d in 1:nd-1, y in 1:nyy], soc_h2[1,d+1,y] == soc_h2[nh+1,d,y] * (1 - h2tank.η_self * controller.Δh) - (p_elyz_E[1,σ_td[d+1,σ_u[y]],y] * elyz.η_E_H2 + p_fc_E[1,σ_td[d+1,σ_u[y]],y] / fc.η_H2_E) * controller.Δh
    # TES dynamic by td
    [h in 1:nh, d in 1:nd, y in 1:nyy], soc_tes[h+1,d,y] == soc_tes[h,d,y] * (1 - tes.η_self * controller.Δh)  - (p_tes_ch[h,σ_td[d,σ_u[y]],y] * tes.η_ch + p_tes_dch[h,σ_td[d,σ_u[y]],y] / tes.η_dch) * controller.Δh
    # TES coupling between td
    [d in 1:nd-1, y in 1:nyy], soc_tes[1,d+1,y] == soc_tes[nh+1,d,y] * (1 - tes.η_self * controller.Δh) - (p_tes_ch[1,σ_td[d+1,σ_u[y]],y] * tes.η_ch  + p_tes_dch[1,σ_td[d+1,σ_u[y]],y] / tes.η_dch) * controller.Δh
    # # Elyz dynamic by td
    [h in 1:nh, d in 1:nd, y in 1:nyy], soh_elyz[h+1,d,y] == soh_elyz[h,d,y] - designer.Δy * δ_elyz[h,σ_td[d,σ_u[y]],y] * controller.Δh / elyz.nHoursMax
    # Elyz coupling between td
    [d in 1:nd-1, y in 1:nyy], soh_elyz[1,d+1,y] == soh_elyz[nh+1,d,y] - designer.Δy * δ_elyz[1,σ_td[d+1,σ_u[y]],y] * controller.Δh / elyz.nHoursMax
    # FC dynamic by td
    [h in 1:nh, d in 1:nd, y in 1:nyy], soh_fc[h+1,d,y] == soh_fc[h,d,y] - designer.Δy * δ_fc[h,σ_td[d,σ_u[y]],y] * controller.Δh / fc.nHoursMax
    # FC coupling between td
    [d in 1:nd-1, y in 1:nyy], soh_fc[1,d+1,y] == soh_fc[nh+1,d,y] - designer.Δy * δ_fc[1,σ_td[d+1,σ_u[y]],y] * controller.Δh / fc.nHoursMax
    # Recourse
    [h in 1:nh, td in 1:ntd, y in 1:nyy], ld_E_bytd[h,td,σ_u[y]] - pMax_pv[y] * pv_bytd[h,td,σ_u[y]] <= p_g_out[h,td,y] + p_g_in[h,td,y] + p_liion_ch[h,td,y] + p_liion_dch[h,td,y] + p_elyz_E[h,td,y] + p_fc_E[h,td,y] + p_heater_E[h,td,y]
    [h in 1:nh, td in 1:ntd, y in 1:nyy], ld_H_bytd[h,td,σ_u[y]] <= - elyz.η_E_H * p_elyz_E[h,td,y] +  fc.η_H2_H / fc.η_H2_E * p_fc_E[h,td,y] -  heater.η_E_H * p_heater_E[h,td,y] + p_tes_ch[h,td,y] + p_tes_dch[h,td,y]
    end)
    # Investment constraints dynamics
    @constraints(m, begin
    # Dynamics
    # PV
    [y in 1:nyy], pMax_pv[y+1] - pMax_pv[y] <= M1 * δ_inv_pv[y]
    [y in 1:nyy], pMax_pv[y] - pMax_pv[y+1] <= M1 * δ_inv_pv[y]
    [y in 1:nyy], pMax_pv[y+1] - r_pv[y] <= M1 * (1 - δ_inv_pv[y])
    [y in 1:nyy], r_pv[y] - pMax_pv[y+1] <= M1 * (1 - δ_inv_pv[y])
    # Liion
    # E_liion
    # if δ_inv_liion = 0 (r_liion = 0) then E_liion[y+1] =  E_liion[y]
    [y in 1:nyy], E_liion[y+1] - E_liion[y] <= M1 * δ_inv_liion[y]
    [y in 1:nyy], E_liion[y] - E_liion[y+1] <= M1 * δ_inv_liion[y]
    # else E_liion[y+1] = r_liion[y] end
    [y in 1:nyy], E_liion[y+1] - r_liion[y] <= M1 * (1 - δ_inv_liion[y])
    [y in 1:nyy], r_liion[y] - E_liion[y+1] <= M1 * (1 - δ_inv_liion[y])
    # SOC_liion
    [y in 1:nyy-1], soc_liion[1,1,y+1] - soc_liion[nh+1,nd,y] <= M1 * δ_inv_liion[y]
    [y in 1:nyy-1], soc_liion[nh+1,nd,y] - soc_liion[1,1,y+1] <= M1 * δ_inv_liion[y]
    [y in 1:nyy-1], soc_liion[1,1,y+1] - liion.α_soc_max * r_liion[y] <= M1 * (1 - δ_inv_liion[y])
    [y in 1:nyy-1], liion.α_soc_max * r_liion[y] - soc_liion[1,1,y+1] <= M1 * (1 - δ_inv_liion[y])
    # SOH_liion
    [y in 1:nyy-1], soh_liion[1,1,y+1] - soh_liion[nh+1,nd,y] <= M1 * δ_inv_liion[y]
    [y in 1:nyy-1], soh_liion[nh+1,nd,y] - soh_liion[1,1,y+1] <= M1 * δ_inv_liion[y]
    [y in 1:nyy-1], soh_liion[1,1,y+1] - r_liion[y] <= M1 * (1 - δ_inv_liion[y])
    [y in 1:nyy-1], r_liion[y] - soh_liion[1,1,y+1] <= M1 * (1 - δ_inv_liion[y])
    # H2 hub
    # E_h2
    [y in 1:nyy], E_h2[y+1] - E_h2[y] <= M2 * δ_inv_tank[y]
    [y in 1:nyy], E_h2[y] - E_h2[y+1] <= M2 * δ_inv_tank[y]
    [y in 1:nyy], E_h2[y+1] - r_tank[y] <= M2 * (1 - δ_inv_tank[y])
    [y in 1:nyy], r_tank[y] - E_h2[y+1] <= M2 * (1 - δ_inv_tank[y])
    # SOC_h2
    [y in 1:nyy-1], soc_h2[1,1,y+1] - soc_h2[nh+1,nd,y] <= M2 * δ_inv_tank[y]
    [y in 1:nyy-1], soc_h2[nh+1,nd,y] - soc_h2[1,1,y+1] <= M2 * δ_inv_tank[y]
    [y in 1:nyy-1], soc_h2[1,1,y+1] - h2tank.α_soc_max * r_tank[y] <= M2 * (1 - δ_inv_tank[y])
    [y in 1:nyy-1], h2tank.α_soc_max * r_tank[y] - soc_h2[1,1,y+1] <= M2 * (1 - δ_inv_tank[y])
    # pMax_elyz
    [y in 1:nyy], pMax_elyz[y+1] - pMax_elyz[y] <= M3 * δ_inv_elyz[y]
    [y in 1:nyy], pMax_elyz[y] - pMax_elyz[y+1] <= M3 * δ_inv_elyz[y]
    [y in 1:nyy], pMax_elyz[y+1] - r_elyz[y] <= M3 * (1 - δ_inv_elyz[y])
    [y in 1:nyy], r_elyz[y] - pMax_elyz[y+1] <= M3 * (1 - δ_inv_elyz[y])
    # SOH_elyz
    [y in 1:nyy-1], soh_elyz[1,1,y+1] - soh_elyz[nh+1,nd,y] <= M3 * δ_inv_elyz[y]
    [y in 1:nyy-1], soh_elyz[nh+1,nd,y] - soh_elyz[1,1,y+1] <= M3 * δ_inv_elyz[y]
    [y in 1:nyy-1], soh_elyz[1,1,y+1] - 1. <= M3 * (1 - δ_inv_elyz[y])
    [y in 1:nyy-1], 1. - soh_elyz[1,1,y+1] <= M3 * (1 - δ_inv_elyz[y])
    # pMax_fc
    [y in 1:nyy], pMax_fc[y+1] - pMax_fc[y] <= M3 * δ_inv_fc[y]
    [y in 1:nyy], pMax_fc[y] - pMax_fc[y+1] <= M3 * δ_inv_fc[y]
    [y in 1:nyy], pMax_fc[y+1] - r_fc[y] <= M3 * (1 - δ_inv_fc[y])
    [y in 1:nyy], r_fc[y] - pMax_fc[y+1] <= M3 * (1 - δ_inv_fc[y])
    # SOH_fc
    [y in 1:nyy-1], soh_fc[1,1,y+1] - soh_fc[nh+1,nd,y] <= M3 * δ_inv_fc[y]
    [y in 1:nyy-1], soh_fc[nh+1,nd,y] - soh_fc[1,1,y+1] <= M3 * δ_inv_fc[y]
    [y in 1:nyy-1], soh_fc[1,1,y+1] - 1. <= M3 * (1 - δ_inv_fc[y])
    [y in 1:nyy-1], 1. - soh_fc[1,1,y+1] <= M3 * (1 - δ_inv_fc[y])
    # TES
    # E_tes
    [y in 1:nyy], E_tes[y+1] - E_tes[y] <= M1 * δ_inv_tes[y]
    [y in 1:nyy], E_tes[y] - E_tes[y+1] <= M1 * δ_inv_tes[y]
    [y in 1:nyy], E_tes[y+1] - r_tes[y] <= M1 * (1 - δ_inv_tes[y])
    [y in 1:nyy], r_tes[y] - E_tes[y+1] <= M1 * (1 - δ_inv_tes[y])
    # SOC_tes
    [y in 1:nyy-1], soc_tes[1,1,y+1] - soc_tes[nh+1,nd,y] <= M1 * δ_inv_tes[y]
    [y in 1:nyy-1], soc_tes[nh+1,nd,y] - soc_tes[1,1,y+1] <= M1 * δ_inv_tes[y]
    [y in 1:nyy-1], soc_tes[1,1,y+1] - tes.α_soc_max * r_tes[y] <= M1 * (1 - δ_inv_tes[y])
    [y in 1:nyy-1], tes.α_soc_max * r_tes[y] - soc_tes[1,1,y+1] <= M1 * (1 - δ_inv_tes[y])
    end)
    # Initial and final conditions
    @constraints(m, begin
    # PV
    pMax_pv[1] == pv.powerMax[1]
    # Liion
    E_liion[1] == liion.Erated[1]
    soh_liion[1,1,1] == liion.soh[1,1] * liion.Erated[1]
    soc_liion[1,1,1] == liion.soc[1,1] * liion.Erated[1]
    [y in 1:nyy], soc_liion[1,1,y] == soc_liion[nh+1,nd,y]
    # TES
    E_tes[1] == tes.Erated[1]
    soc_tes[1,1,1] == tes.soc[1,1] * tes.Erated[1]
    [y in 1:nyy], soc_tes[1,1,y] == soc_tes[nh+1,nd,y]
    # H2
    E_h2[1] == h2tank.Erated[1]
    soc_h2[1,1,1] == h2tank.soc[1,1] * h2tank.Erated[1]
    [y in 1:nyy], soc_h2[1,1,y] == soc_h2[nh+1,nd,y]
    # Elyz
    pMax_elyz[1] == elyz.powerMax[1]
    soh_elyz[1,1,1] == elyz.soh[1,1]
    # FC
    pMax_fc[1] == fc.powerMax[1]
    soh_fc[1,1,1] == fc.soh[1,1]
    end)
    # Grid constraints
    @constraints(m, begin
    [h in 1:nh, td in 1:ntd, y in 2:nyy], p_g_in[h,td,y] <= (1. - grid.τ_power) * maximum(ld_E_bytd[:,:,σ_u[y]] .+ ld_H_bytd[:,:,σ_u[y]] / heater.η_E_H)
    [y in 2:nyy], sum(p_g_in[h,td,y] for h in 1:nh, td in 1:ntd) <= (1. - grid.τ_energy) * sum(ld_E_bytd[h,td,σ_u[y]] for h in 1:nh, td in 1:ntd)
    end)
    # Objective
    @objective(m, Min, sum( γ[σ_u[y]] * (pv.C_pv[σ_u[y]] * r_pv[y] + liion.C_liion[σ_u[y]] * r_liion[y] + h2tank.C_tank[σ_u[y]] * r_tank[y] +
    elyz.C_elyz[σ_u[y]] * r_elyz[y] + fc.C_fc[σ_u[y]] * r_fc[y] + tes.C_tes[σ_u[y]] * r_tes[y] +
    designer.Δy * sum(n_bytd[td] * (p_g_in[h,td,y] * grid.C_grid_in[h, σ_u[y]] + p_g_out[h,td,y] * grid.C_grid_out[h, σ_u[y]]) * controller.Δh  for h=1:nh, td in 1:ntd) ) for y=1:nyy))
    # Solve
    optimize!(m)
    # Outputs
    # Operation controls
    controller.u = (
    u_liion =  value.(p_liion_ch + p_liion_dch),
    u_elyz = value.(p_elyz_E),
    u_fc= value.(p_fc_E),
    u_tes = value.(p_tes_ch + p_tes_dch),
    u_heater = value.(p_heater_E),
    )
    # Investment controls
    designer.u = (
    u_pv = value.(r_pv),
    u_liion = value.(r_liion),
    u_tank = value.(r_tank),
    u_elyz = value.(r_elyz),
    u_fc = value.(r_fc),
    u_tes = value.(r_tes),
    )
    # Operation states
    controller.x = (
    soc_liion = value.(soc_liion),
    soc_h2 = value.(soc_h2),
    soh_elyz = value.(soh_elyz),
    soh_fc = value.(soh_fc),
    soc_tes = value.(soc_tes),
    soh_liion = value.(soh_liion),
    )
    # Investment states
    designer.x = (
    pMax_pv = value.(pMax_pv),
    E_liion = value.(E_liion),
    E_h2 = value.(E_h2),
    pMax_elyz = value.(pMax_elyz),
    pMax_fc = value.(pMax_fc),
    E_tes = value.(E_tes),
    )
    # Log
    controller.log = designer.log = termination_status(m)
    # Model
    controller.model = designer.model = m
end

#                                   Offline functions
#_______________________________________________________________________________

#                                   Online functions
#_______________________________________________________________________________
