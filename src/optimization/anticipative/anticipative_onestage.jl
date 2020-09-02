#=
    Anticipative one-stage designer with anticipative controller
=#

mutable struct AnticipativeOneStageDesigner <: AbstractDesigner
    u::NamedTuple
    horizon::Int64
    model::JuMP.Model
    AnticipativeOneStageDesigner() = new()
end

#### Models ####
# Simple
function onestage_milp_model(ld::Load, pv::Source, liion::Liion,
  controller::AnticipativeController, designer::AnticipativeOneStageDesigner,
  grid::Grid, ω_optim::Scenarios, parameters::NamedTuple)

    # Parameters
    crf_liion = (parameters.τ * (parameters.τ + 1.) ^ liion.lifetime) / ((parameters.τ + 1.) ^ liion.lifetime - 1.)
    crf_pv = (parameters.τ * (parameters.τ + 1.) ^ pv.lifetime) / ((parameters.τ + 1.) ^ pv.lifetime - 1.)

    # Sets
    nh = length(1:parameters.Δh:controller.horizon) # Number of hours
    ns = size(ω_optim.ld_E,2) # Number of one-year scenarios

    # Model definition
    m = Model(CPLEX.Optimizer)

    # Variables
    # Operation control variables
    @variables(m, begin
    # Liion
     p_liion_ch[1:nh,1:ns] <= 0.
     p_liion_dch[1:nh,1:ns] >= 0.
     # Recourse
     p_g_out[1:nh,1:ns] <= 0.
     p_g_in[1:nh,1:ns] >= 0.
    end)

    # Investment control variables
    @variables(m, begin
    r_liion >= 0.
    r_pv >= 0.
    end)

    # Operation state variables
    @variables(m, begin
    soc_liion[1:nh+1,1:ns]
    end)

    # Operation constraints bounds
    @constraints(m, begin
    # Controls
    # Liion
    [h in 1:nh, s in 1:ns], p_liion_dch[h,s] <= liion.α_p_dch * r_liion
    [h in 1:nh, s in 1:ns], p_liion_ch[h,s] >= -liion.α_p_ch * r_liion

    # State
    # Liion
    [h in 1:nh+1, s in 1:ns], soc_liion[h,s] <= liion.α_soc_max * r_liion
    [h in 1:nh+1, s in 1:ns], soc_liion[h,s] >= liion.α_soc_min * r_liion
    end)

    # Investment constraints bounds
    @constraints(m, begin
    # Controls
    r_liion <= 1000.
    r_pv <= 1000.
    end)

    # Operation constraints dynamics and recourse
    @constraints(m, begin
    # Dynamics
    # Liion
    [h in 1:nh, s in 1:ns], soc_liion[h+1,s] == soc_liion[h,s] * (1 - liion.η_self * parameters.Δh) - (p_liion_ch[h,s] * liion.η_ch + p_liion_dch[h,s] / liion.η_dch) * parameters.Δh

    # Recourse
    [h in 1:nh, s in 1:ns], ω_optim.ld_E[h,s] - r_pv * ω_optim.pv_E[h,s] <= p_g_out[h,s] + p_g_in[h,s] + p_liion_ch[h,s] + p_liion_dch[h,s]
    end)

    # Initial and final conditions
    @constraints(m, begin
    [s in 1:ns], soc_liion[1,s] == liion.soc[1] * r_liion
    [s in 1:ns], soc_liion[end,s] == soc_liion[1,s]
    end)

    # Grid constraints
    @constraints(m, begin
    power_constraint, p_g_in .<= (1. - grid.τ_power) * maximum(ω_optim.ld_E)
    self_constraint[s in 1:ns], sum(p_g_in[h,s] for h in 1:nh) <= (1. - grid.τ_energy) * sum(ω_optim.ld_E[h,s] for h in 1:nh)
    end)

    # CAPEX
    capex = @expression(m, crf_pv * ω_optim.C_pv[1] * r_pv + crf_liion * ω_optim.C_liion[1] * r_liion)

    # OPEX
    opex = @expression(m, sum((p_g_in[h,s] * ω_optim.C_grid_in[h,s] + p_g_out[h,s] * ω_optim.C_grid_out[h,s]) * parameters.Δh  for h in 1:nh, s in 1:ns) / ns)

    # Objective
    @objective(m, Min, capex + opex)

    return m
end
# Multi-energy
function onestage_milp_model(ld::Load, pv::Source, liion::Liion, h2tank::H2Tank,
   elyz::Electrolyzer, fc::FuelCell, tes::ThermalSto, heater::Heater,
    controller::AnticipativeController, designer::AnticipativeOneStageDesigner,
    grid::Grid, ω_optim::Scenarios, parameters::NamedTuple)
    # Parameters
    crf_pv = (parameters.τ * (parameters.τ + 1.) ^ pv.lifetime) / ((parameters.τ + 1.) ^ pv.lifetime - 1.)
    crf_liion = (parameters.τ * (parameters.τ + 1.) ^ liion.lifetime) / ((parameters.τ + 1.) ^ liion.lifetime - 1.)
    crf_tank = (parameters.τ * (parameters.τ + 1.) ^ h2tank.lifetime) / ((parameters.τ + 1.) ^ h2tank.lifetime - 1.)
    crf_elyz = (parameters.τ * (parameters.τ + 1.) ^ elyz.lifetime) / ((parameters.τ + 1.) ^ elyz.lifetime - 1.)
    crf_fc = (parameters.τ * (parameters.τ + 1.) ^ fc.lifetime) / ((parameters.τ + 1.) ^ fc.lifetime - 1.)
    crf_tes = (parameters.τ * (parameters.τ + 1.) ^ tes.lifetime) / ((parameters.τ + 1.) ^ tes.lifetime - 1.)
    # Sets
    nhh = length(1:parameters.Δh:controller.horizon) # Number of hours
    ns = size(ω_optim.ld_E,2) # Number of one-year scenarios
    # Model definition
    m = Model(CPLEX.Optimizer)
    # Variables
    # Operation decisions variables
    @variables(m, begin
    # Liion
    p_liion_ch[1:nhh, 1:ns] <= 0.
    p_liion_dch[1:nhh, 1:ns] >= 0.
    # H2 hub
    p_elyz_E[1:nhh, 1:ns] <= 0.
    p_fc_E[1:nhh, 1:ns] >= 0.
    # TES
    p_tes_ch[1:nhh, 1:ns] <= 0.
    p_tes_dch[1:nhh, 1:ns] >= 0.
    # Heater
    p_heater_E[1:nhh, 1:ns] <= 0.
    # Recourse
    p_g_out[1:nhh, 1:ns] <= 0.
    p_g_in[1:nhh, 1:ns] >= 0.
    end)
    # Operation state variables
    @variables(m, begin
    soc_liion[1:nhh+1, 1:ns]
    soc_h2[1:nhh+1, 1:ns]
    soc_tes[1:nhh+1, 1:ns]
    end)
    # Investment decisions variables
    @variables(m, begin
    0 <= r_pv <= 1000
    0 <= r_liion <= 1000
    0 <= r_tank <= 50000
    0 <= r_elyz <= 50
    0 <= r_fc <= 50
    0 <= r_tes <= 1000
    end)
    # Operation constraints bounds
    @constraints(m, begin
    # Controls
    # Liion
    [h in 1:nhh, s in 1:ns], p_liion_dch[h,s] <= liion.α_p_dch * r_liion
    [h in 1:nhh, s in 1:ns], p_liion_ch[h,s] >= -liion.α_p_ch * r_liion
    # TES
    [h in 1:nhh, s in 1:ns], p_tes_dch[h,s] <= tes.α_p_dch * r_tes
    [h in 1:nhh, s in 1:ns], p_tes_ch[h,s] >= -tes.α_p_ch * r_tes
    # H2 hub
    [h in 1:nhh, s in 1:ns], p_fc_E[h,s] <= r_fc
    [h in 1:nhh, s in 1:ns], p_elyz_E[h,s] >= -r_elyz
    [h in 1:nhh, s in 1:ns], p_heater_E[h,s] >= -heater.powerMax[1]
    # State
    # Liion
    [h in 1:nhh+1, s in 1:ns], soc_liion[h,s] <= liion.α_soc_max * r_liion
    [h in 1:nhh+1, s in 1:ns], soc_liion[h,s] >= liion.α_soc_min * r_liion
    # H2 hub
    [h in 1:nhh+1, s in 1:ns], soc_h2[h,s] <= h2tank.α_soc_max * r_tank
    [h in 1:nhh+1, s in 1:ns], soc_h2[h,s] >= h2tank.α_soc_min * r_tank
    # TES
    [h in 1:nhh+1, s in 1:ns], soc_tes[h,s] <= tes.α_soc_max * r_tes
    [h in 1:nhh+1, s in 1:ns], soc_tes[h,s] >= tes.α_soc_min * r_tes
    end)
    # Operation constraints dynamics and recourse
    @constraints(m, begin
    # Dynamics
    # Liion
    [h in 1:nhh, s in 1:ns], soc_liion[h+1,s] == soc_liion[h,s] * (1. - liion.η_self * parameters.Δh) - (p_liion_ch[h,s] * liion.η_ch + p_liion_dch[h,s] / liion.η_dch) * parameters.Δh
    # H2 hub
    [h in 1:nhh, s in 1:ns], soc_h2[h+1,s] == soc_h2[h,s] * (1. - h2tank.η_self * parameters.Δh) - (p_elyz_E[h,s] * elyz.η_E_H2 + p_fc_E[h,s] / fc.η_H2_E) * parameters.Δh
    # TES
    [h in 1:nhh, s in 1:ns], soc_tes[h+1,s] == soc_tes[h,s] * (1. - tes.η_self * parameters.Δh) - (p_tes_ch[h,s] * tes.η_ch + p_tes_dch[h,s] / tes.η_dch) * parameters.Δh
    # Recourse
    [h in 1:nhh, s in 1:ns], ω_optim.ld_E[h,s] - r_pv * ω_optim.pv_E[h,s] <= p_g_out[h,s] + p_g_in[h,s] + p_liion_ch[h,s] + p_liion_dch[h,s] + p_elyz_E[h,s] + p_fc_E[h,s] + p_heater_E[h,s]
    [h in 1:nhh, s in 1:ns], ω_optim.ld_H[h,s] <= - elyz.η_E_H * p_elyz_E[h,s] +  fc.η_H2_H / fc.η_H2_E * p_fc_E[h,s] -  heater.η_E_H * p_heater_E[h,s] + p_tes_ch[h,s] + p_tes_dch[h,s]
    end)
    # Initial and final conditions
    @constraints(m, begin
    # Liion
    [s in 1:ns], soc_liion[1,s] == liion.soc[1,1] * r_liion
    [s in 1:ns], soc_liion[end,s] >= soc_liion[1,s]
    # H2
    [s in 1:ns], soc_h2[1,s] == h2tank.soc[1,1] * r_tank
    [s in 1:ns], soc_h2[end,s] >= soc_h2[1,s]
    # TES
    [s in 1:ns], soc_tes[1,s] == tes.soc[1,1] * r_tes
    [s in 1:ns], soc_tes[end,s] >= soc_tes[1,s]
    end)
    # Grid constraints
    @constraints(m, begin
    power_constraint, p_g_in .<= (1. - grid.τ_power) * maximum(ω_optim.ld_E .+ ω_optim.ld_H ./ heater.η_E_H) # seems to be faster in vectorized form
    self_constraint[s in 1:ns], sum(p_g_in[h,s] for h in 1:nhh) <= (1. - grid.τ_energy) * sum(ω_optim.ld_E[h,s] .+ ω_optim.ld_H[h,s] ./ heater.η_E_H for h in 1:nhh)
    end)
    # CAPEX
    capex = @expression(m, crf_pv * ω_optim.C_pv[1] * r_pv +
    crf_liion * ω_optim.C_liion[1] * r_liion +
    crf_tank * ω_optim.C_tank[1] * r_tank +
    crf_elyz * ω_optim.C_elyz[1] * r_elyz +
    crf_fc * ω_optim.C_fc[1] * r_fc +
    crf_tes * ω_optim.C_tes[1] * r_tes)
    # OPEX
    opex = @expression(m, sum((p_g_in[h,s] * ω_optim.C_grid_in[h,s] + p_g_out[h,s] * ω_optim.C_grid_out[h,s]) * parameters.Δh  for h in 1:nhh, s in 1:ns) / ns)
    # Objective
    @objective(m, Min, capex + opex)
    return m
end
# Multi-energy with clustering
function onestage_milp_model(ld::Load, pv::Source, liion::Liion, h2tank::H2Tank,
   elyz::Electrolyzer, fc::FuelCell, tes::ThermalSto, heater::Heater,
    controller::AnticipativeController, designer::AnticipativeOneStageDesigner,
    grid::Grid, ω_optim::Scenarios, ω_td::ClusteredScenarios, parameters::NamedTuple)
    # Parameters
    crf_pv = (parameters.τ * (parameters.τ + 1.) ^ pv.lifetime) / ((parameters.τ + 1.) ^ pv.lifetime - 1.)
    crf_liion = (parameters.τ * (parameters.τ + 1.) ^ liion.lifetime) / ((parameters.τ + 1.) ^ liion.lifetime - 1.)
    crf_tank = (parameters.τ * (parameters.τ + 1.) ^ h2tank.lifetime) / ((parameters.τ + 1.) ^ h2tank.lifetime - 1.)
    crf_elyz = (parameters.τ * (parameters.τ + 1.) ^ elyz.lifetime) / ((parameters.τ + 1.) ^ elyz.lifetime - 1.)
    crf_fc = (parameters.τ * (parameters.τ + 1.) ^ fc.lifetime) / ((parameters.τ + 1.) ^ fc.lifetime - 1.)
    crf_tes = (parameters.τ * (parameters.τ + 1.) ^ tes.lifetime) / ((parameters.τ + 1.) ^ tes.lifetime - 1.)
    # Sets
    nhh = size(ω_td.clusters.ld_E,1) # number of hours by td (=24)
    ntd = size(ω_td.clusters.ld_E,2) # number of td
    nd = size(ω_td.σ,1) # number of days over the year (length of the sequence)
    ns = size(ω_td.clusters.ld_E,3)  # Number of one-year scenarios
    # Model definition
    m = Model(CPLEX.Optimizer)
    # Variables
    # Operation control variables
    @variables(m, begin
    # Liion
     p_liion_ch[1:nhh, 1:ntd, 1:ns] <= 0.
     p_liion_dch[1:nhh, 1:ntd, 1:ns] >= 0.
     # H2 hub
     p_elyz_E[1:nhh, 1:ntd, 1:ns] <= 0.
     p_fc_E[1:nhh, 1:ntd, 1:ns] >= 0.
     # TES
     p_tes_ch[1:nhh, 1:ntd, 1:ns] <= 0.
     p_tes_dch[1:nhh, 1:ntd, 1:ns] >= 0.
     # Heater
     p_heater_E[1:nhh, 1:ntd, 1:ns] <= 0.
     # Recourse
     p_g_out[1:nhh, 1:ntd, 1:ns] <= 0.
     p_g_in[1:nhh, 1:ntd, 1:ns] >= 0.
    end)
    # Operation state variables
    @variables(m, begin
    soc_liion[1:nhh+1, 1:nd, 1:ns]
    soc_h2[1:nhh+1, 1:nd, 1:ns]
    soc_tes[1:nhh+1, 1:nd, 1:ns]
    end)
    # Investment decisions variables
    @variables(m, begin
    0 <= r_pv <= 1000
    0 <= r_liion <= 1000
    0 <= r_tank <= 50000
    0 <= r_elyz <= 50
    0 <= r_fc <= 50
    0 <= r_tes <= 1000
    end)
    # Operation constraints bounds
    @constraints(m, begin
    # Controls
    # Liion
    [h in 1:nhh, td in 1:ntd, s in 1:ns], p_liion_dch[h,td,s] <= liion.α_p_dch * r_liion
    [h in 1:nhh, td in 1:ntd, s in 1:ns], p_liion_ch[h,td,s] >= -liion.α_p_ch * r_liion
    # TES
    [h in 1:nhh, td in 1:ntd, s in 1:ns], p_tes_dch[h,td,s] <= tes.α_p_dch * r_tes
    [h in 1:nhh, td in 1:ntd, s in 1:ns], p_tes_ch[h,td,s] >= -tes.α_p_ch * r_tes
    # H2 hub
    [h in 1:nhh, td in 1:ntd, s in 1:ns], p_fc_E[h,td,s] <= r_fc
    [h in 1:nhh, td in 1:ntd, s in 1:ns], p_elyz_E[h,td,s] >= -r_elyz
    [h in 1:nhh, td in 1:ntd, s in 1:ns], p_heater_E[h,td,s] >= -heater.powerMax[1]
    # State
    # Liion
    [h in 1:nhh+1, d in 1:nd, s in 1:ns], soc_liion[h,d,s] <= liion.α_soc_max * r_liion
    [h in 1:nhh+1, d in 1:nd, s in 1:ns], soc_liion[h,d,s] >= liion.α_soc_min * r_liion
    # H2 hub
    [h in 1:nhh+1, d in 1:nd, s in 1:ns], soc_h2[h,d,s] <= h2tank.α_soc_max * r_tank
    [h in 1:nhh+1, d in 1:nd, s in 1:ns], soc_h2[h,d,s] >= h2tank.α_soc_min * r_tank
    # TES
    [h in 1:nhh+1, d in 1:nd, s in 1:ns], soc_tes[h,d,s] <= tes.α_soc_max * r_tes
    [h in 1:nhh+1, d in 1:nd, s in 1:ns], soc_tes[h,d,s] >= tes.α_soc_min * r_tes
    end)
    # Operation constraints dynamics and recourse
    @constraints(m, begin
    # Dynamics
    # Liion
    [h in 1:nhh, d in 1:nd, s in 1:ns], soc_liion[h+1,d,s] == soc_liion[h,d,s] * (1. - liion.η_self * parameters.Δh) - (p_liion_ch[h,ω_td.σ[d,s],s] * liion.η_ch + p_liion_dch[h,ω_td.σ[d,s],s] / liion.η_dch) * parameters.Δh
    [d in 1:nd-1, s in 1:ns], soc_liion[1,d+1,s] == soc_liion[nhh+1,d,s] * (1. - liion.η_self * parameters.Δh) - (p_liion_ch[1,ω_td.σ[d+1,s],s] * liion.η_ch + p_liion_dch[1,ω_td.σ[d+1,s],s] / liion.η_dch) * parameters.Δh
    # H2 hub
    [h in 1:nhh, d in 1:nd, s in 1:ns], soc_h2[h+1,d,s] == soc_h2[h,d,s] * (1. - h2tank.η_self * parameters.Δh) - (p_elyz_E[h,ω_td.σ[d,s],s] * elyz.η_E_H2 + p_fc_E[h,ω_td.σ[d,s],s] / fc.η_H2_E) * parameters.Δh
    [d in 1:nd-1, s in 1:ns], soc_h2[1,d+1,s] == soc_h2[nhh+1,d,s] * (1. - h2tank.η_self * parameters.Δh) - (p_elyz_E[1,ω_td.σ[d+1,s],s] * elyz.η_E_H2 + p_fc_E[1,ω_td.σ[d+1,s],s] / fc.η_H2_E) * parameters.Δh
    # TES
    [h in 1:nhh, d in 1:nd, s in 1:ns], soc_tes[h+1,d,s] == soc_tes[h,d,s] * (1. - tes.η_self * parameters.Δh) - (p_tes_ch[h,ω_td.σ[d,s],s] * tes.η_ch + p_tes_dch[h,ω_td.σ[d,s],s] / tes.η_dch) * parameters.Δh
    [d in 1:nd-1, s in 1:ns], soc_tes[1,d+1,s] == soc_tes[nhh+1,d,s] * (1. - tes.η_self * parameters.Δh) - (p_tes_ch[1,ω_td.σ[d+1,s],s] * tes.η_ch + p_tes_dch[1,ω_td.σ[d+1,s],s] / tes.η_dch) * parameters.Δh
    # Recourse
    [h in 1:nhh, td in 1:ntd, s in 1:ns], ω_td.clusters.ld_E[h,td,s] - r_pv * ω_td.clusters.pv_E[h,td,s] <= p_g_out[h,td,s] + p_g_in[h,td,s] + p_liion_ch[h,td,s] + p_liion_dch[h,td,s] + p_elyz_E[h,td,s] + p_fc_E[h,td,s] + p_heater_E[h,td,s]
    [h in 1:nhh, td in 1:ntd, s in 1:ns], ω_td.clusters.ld_H[h,td,s] <= - elyz.η_E_H * p_elyz_E[h,td,s] +  fc.η_H2_H / fc.η_H2_E * p_fc_E[h,td,s] -  heater.η_E_H * p_heater_E[h,td,s] + p_tes_ch[h,td,s] + p_tes_dch[h,td,s]
    end)
    # Initial and final conditions
    @constraints(m, begin
    # Liion
    [s in 1:ns], soc_liion[1,1,s] == liion.soc[1,1] * r_liion
    [s in 1:ns], soc_liion[end,end,s] >= soc_liion[1,1,s]
    # H2
    [s in 1:ns], soc_h2[1,1,s] == h2tank.soc[1,1] * r_tank
    [s in 1:ns], soc_h2[end,end,s] >= soc_h2[1,1,s]
    # TES
    [s in 1:ns], soc_tes[1,1,s] == tes.soc[1,1] * r_tes
    [s in 1:ns], soc_tes[end,end,s] >= soc_tes[1,1,s]
    end)
    # Grid constraints
    @constraints(m, begin
    power_constraint, p_g_in[1:nhh,1:ntd,1:ns] .<= (1. - grid.τ_power) * maximum(ω_optim.ld_E + ω_optim.ld_H / heater.η_E_H)
    self_constraint[s in 1:ns], sum(ω_td.nby[td,s] * sum(p_g_in[h,td,s] for h in 1:nhh) for td in 1:ntd) <= (1. - grid.τ_energy) * sum(ω_optim.ld_E[h,s] + ω_optim.ld_H[h,s] / heater.η_E_H for h in 1:nhh)
    end)
    # CAPEX
    capex = @expression(m, crf_pv * ω_optim.C_pv[1] * r_pv +
    crf_liion * ω_optim.C_liion[1] * r_liion +
    crf_tank * ω_optim.C_tank[1] * r_tank +
    crf_elyz * ω_optim.C_elyz[1] * r_elyz +
    crf_fc * ω_optim.C_fc[1] * r_fc +
    crf_tes * ω_optim.C_tes[1] * r_tes)
    # OPEX
    opex = @expression(m, sum(sum(ω_td.nby[td,s] * sum((p_g_in[h,td,s] * ω_optim.C_grid_in[h,s] + p_g_out[h,td,s] * ω_optim.C_grid_out[h,s]) * parameters.Δh  for h in 1:nhh) for td in 1:ntd) for s in 1:ns) / ns)
    # Objective
    @objective(m, Min, capex + opex)
    return m
end

#### Offline functions ####
# Simple
function offline_optimization(ld::Load, pv::Source, liion::Liion,
     controller::AnticipativeController, designer::AnticipativeOneStageDesigner,
     grid::Grid, ω_optim::Scenarios, parameters::NamedTuple)
     # Parameters
     ny = size(ld.power_E,2) # number of simulation years
     ns = size(ld.power_E,3) # number of scenarios

     # Selection of one-year scenarios from the optimization dataset
     ω_milp = Scenarios(ω_optim.timestamp,ω_optim.ld_E[:,:,1],nothing,ω_optim.pv_E[:,:,1],
     ω_optim.C_pv[:,1],ω_optim.C_liion[:,1],nothing,nothing,nothing,nothing,nothing,
     ω_optim.C_grid_in[:,:,1],ω_optim.C_grid_out[:,:,1])

     # Initialize model
     designer.model = controller.model = onestage_milp_model(
     ld, pv, liion, controller, designer, grid, ω_milp, parameters)

     # Compute both investment and operation decisions
     optimize!(designer.model)

     # Formatting variables to simulation
     # Operation decisions
     controller.u = (
     u_liion =  repeat(value.(controller.model[:p_liion_ch] + controller.model[:p_liion_dch]), 1, 1, ns),
     )
     # Investment controls
     designer.u = (
     u_liion = repeat(vcat(value.(designer.model[:r_liion]), zeros(ny-1,1)),1, ns),
     u_pv = repeat(vcat(value.(designer.model[:r_pv]), zeros(ny-1,1)),1, ns),
     )
end
# Multi-energy
function offline_optimization(ld::Load, pv::Source, liion::Liion, h2tank::H2Tank,
   elyz::Electrolyzer, fc::FuelCell, tes::ThermalSto, heater::Heater,
   controller::AnticipativeController, designer::AnticipativeOneStageDesigner,
   grid::Grid, ω_optim::Scenarios, parameters::NamedTuple)
     # Parameters
     ny = size(ld.power_E,2) # number of simulation years
     ns = size(ld.power_E,3) # number of scenarios

     # Selection of one-year scenarios from the optimization dataset
     ω_milp = Scenarios(ω_optim.timestamp, ω_optim.ld_E[:,:,1], ω_optim.ld_H[:,:,1],
     ω_optim.pv_E[:,:,1], ω_optim.C_pv[:,1], ω_optim.C_liion[:,1], ω_optim.C_tes[:,1],
     ω_optim.C_tank[:,1], ω_optim.C_elyz[:,1], ω_optim.C_fc[:,1], ω_optim.C_heater[:,1],
     ω_optim.C_grid_in[:,:,1], ω_optim.C_grid_out[:,:,1])

     # Initialize model
     designer.model = controller.model = onestage_milp_model(
     ld, pv, liion, h2tank, elyz, fc, tes, heater, controller, designer, grid, ω_milp, parameters)

     # Compute both investment and operation decisions
     optimize!(designer.model)

     # Formatting variables to simulation
     # Operation decisions
     controller.u = (
     u_liion =  repeat(value.(controller.model[:p_liion_ch] + controller.model[:p_liion_dch]), 1, 1, ns),
     u_elyz = repeat(value.(controller.model[:p_elyz_E]), 1, 1, ns),
     u_fc= repeat(value.(controller.model[:p_fc_E]), 1, 1, ns),
     u_tes = repeat(value.(controller.model[:p_tes_ch] + controller.model[:p_tes_dch]), 1, 1, ns),
     u_heater = repeat(value.(controller.model[:p_heater_E]), 1, 1, ns),
     )
     # Investment controls
     designer.u = (
     u_pv = repeat(vcat( value.(designer.model[:r_pv]), zeros(ny-1,1)), 1, ns),
     u_liion = repeat(vcat(value.(designer.model[:r_liion]), zeros(ny-1,1)), 1, ns),
     u_tank = repeat(vcat(value.(designer.model[:r_tank]), zeros(ny-1,1)), 1, ns),
     u_elyz = repeat(vcat(value.(designer.model[:r_elyz]), zeros(ny-1,1)), 1, ns),
     u_fc = repeat(vcat(value.(designer.model[:r_fc]), zeros(ny-1,1)), 1, ns),
     u_tes = repeat(vcat(value.(designer.model[:r_tes]), zeros(ny-1,1)), 1, ns),
     )
end
# Multi-energy with clustering
function offline_optimization(ld::Load, pv::Source, liion::Liion, h2tank::H2Tank,
   elyz::Electrolyzer, fc::FuelCell, tes::ThermalSto, heater::Heater,
   controller::AnticipativeController, designer::AnticipativeOneStageDesigner,
   grid::Grid, ω_optim::Scenarios, ntd::Int64, parameters::NamedTuple)
     # Parameters
     ny = size(ld.power_E,2) # number of simulation years
     ns = size(ld.power_E,3) # number of scenarios

     # Selection of one-year scenarios from the optimization dataset
     ω_milp = Scenarios(ω_optim.timestamp, ω_optim.ld_E[:,:,1],ω_optim.ld_H[:,:,1],ω_optim.pv_E[:,:,1],
     ω_optim.C_pv[:,1],ω_optim.C_liion[:,1],ω_optim.C_tes[:,1],ω_optim.C_tank[:,1],
     ω_optim.C_elyz[:,1],ω_optim.C_fc[:,1],ω_optim.C_heater[:,1],
     ω_optim.C_grid_in[:,:,1],ω_optim.C_grid_out[:,:,1])

     # Typical days clustering
     ω_td = clustering_typical_day(ω_milp, ntd)

     # Initialize model
     designer.model = onestage_milp_model(
     ld, pv, liion, h2tank, elyz, fc, tes, heater, controller, designer, grid, ω_milp, ω_td, parameters)

     # Compute both investment and operation decisions
     optimize!(designer.model)

     # Initialize model for the operation
     controller.model = onestage_milp_model(
     ld, pv, liion, h2tank, elyz, fc, tes, heater, controller, designer, grid, ω_milp, parameters)

     # Delete self-sufficiency constraint for the operation
     delete.(controller.model, controller.model[:self_constraint])

     # Fix sizes from investment optimization
     fix.(controller.model[:r_pv], value.(designer.model[:r_pv]), force = true)
     fix.(controller.model[:r_liion], value.(designer.model[:r_liion]), force = true)
     fix.(controller.model[:r_tank], value.(designer.model[:r_tank]), force = true)
     fix.(controller.model[:r_elyz], value.(designer.model[:r_elyz]), force = true)
     fix.(controller.model[:r_fc], value.(designer.model[:r_fc]), force = true)
     fix.(controller.model[:r_tes], value.(designer.model[:r_tes]), force = true)

     # Compute operation decisions with initial scenarios
     optimize!(controller.model)

     # Formatting variables to simulation
     # Operation decisions
     controller.u = (
     u_liion =  repeat(value.(controller.model[:p_liion_ch] + controller.model[:p_liion_dch]), 1, 1, ns),
     u_elyz = repeat(value.(controller.model[:p_elyz_E]), 1, 1, ns),
     u_fc= repeat(value.(controller.model[:p_fc_E]), 1, 1, ns),
     u_tes = repeat(value.(controller.model[:p_tes_ch] + controller.model[:p_tes_dch]), 1, 1, ns),
     u_heater = repeat(value.(controller.model[:p_heater_E]), 1, 1, ns),
     )
     # Investment controls
     designer.u = (
     u_pv = repeat(vcat( value.(designer.model[:r_pv]), zeros(ny-1,1)), 1, ns),
     u_liion = repeat(vcat(value.(designer.model[:r_liion]), zeros(ny-1,1)), 1, ns),
     u_tank = repeat(vcat(value.(designer.model[:r_tank]), zeros(ny-1,1)), 1, ns),
     u_elyz = repeat(vcat(value.(designer.model[:r_elyz]), zeros(ny-1,1)), 1, ns),
     u_fc = repeat(vcat(value.(designer.model[:r_fc]), zeros(ny-1,1)), 1, ns),
     u_tes = repeat(vcat(value.(designer.model[:r_tes]), zeros(ny-1,1)), 1, ns),
     )
end

#### Online functions ####
# Simple
function compute_investment_decisions(y::Int64, s::Int64, ld::Load, pv::Source,
    liion::Liion, grid::Grid, controller::AnticipativeController,
    designer::AnticipativeOneStageDesigner, ω_optim::Scenarios, parameters::NamedTuple)
    ϵ = 0.1
    if liion.soh[end,y,s] < ϵ
        designer.u.u_liion[y,s] = liion.Erated[y,s]
    end
end
# Multi-energy
function compute_investment_decisions(y::Int64, s::Int64, ld::Load, pv::Source,
    liion::Liion, h2tank::H2Tank, elyz::Electrolyzer, fc::FuelCell, tes::ThermalSto,
    heater::Heater, designer::AnticipativeOneStageDesigner, ω_optim::Scenarios, parameters::NamedTuple)
    ϵ = 0.1

    # Liion
    if liion.soh[end,y,s] < ϵ
        designer.u.u_liion[y,s] = liion.Erated[y,s]
    end

    # Electrolyzer
    if elyz.soh[end,y,s] < ϵ
        designer.u.u_elyz[y,s] = elyz.powerMax[y,s]
    end

    # FuelCell
    if fc.soh[end,y,s] < ϵ
        designer.u.u_fc[y,s] = fc.powerMax[y,s]
    end
end
