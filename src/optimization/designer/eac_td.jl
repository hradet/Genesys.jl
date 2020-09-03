#=
    Designer based on the equivalent annual cost (EAC)
=#

mutable struct TypicalDayEACDesigner <: AbstractDesigner
    u::NamedTuple
    horizon::Int64
    ntd::Int64
    model::JuMP.Model
    TypicalDayEACDesigner() = new()
end

#### Models ####
# Multi-energy
function eac_milp_model(ld::Load, pv::Source, liion::Liion, h2tank::H2Tank,
   elyz::Electrolyzer, fc::FuelCell, tes::ThermalSto, heater::Heater,
    designer::TypicalDayEACDesigner, grid::Grid, ω_optim::ClusteredScenarios,
    parameters::NamedTuple)
    # Parameters
    Γ_pv = (parameters.τ * (parameters.τ + 1.) ^ pv.lifetime) / ((parameters.τ + 1.) ^ pv.lifetime - 1.)
    Γ_liion = (parameters.τ * (parameters.τ + 1.) ^ liion.lifetime) / ((parameters.τ + 1.) ^ liion.lifetime - 1.)
    Γ_tank = (parameters.τ * (parameters.τ + 1.) ^ h2tank.lifetime) / ((parameters.τ + 1.) ^ h2tank.lifetime - 1.)
    Γ_elyz = (parameters.τ * (parameters.τ + 1.) ^ elyz.lifetime) / ((parameters.τ + 1.) ^ elyz.lifetime - 1.)
    Γ_fc = (parameters.τ * (parameters.τ + 1.) ^ fc.lifetime) / ((parameters.τ + 1.) ^ fc.lifetime - 1.)
    Γ_tes = (parameters.τ * (parameters.τ + 1.) ^ tes.lifetime) / ((parameters.τ + 1.) ^ tes.lifetime - 1.)

    # Sets
    nh = size(ω_optim.clusters.ld_E,1) # number of hours by td (=24)
    ntd = size(ω_optim.clusters.ld_E,2) # number of td
    nd = size(ω_optim.σ,1) # number of days over the year (length of the sequence)

    # Model definition
    m = Model(CPLEX.Optimizer)

    # Variables
    @variables(m, begin
     # Operation decision variables
     p_liion_ch[1:nh, 1:ntd] <= 0.
     p_liion_dch[1:nh, 1:ntd] >= 0.
     p_elyz_E[1:nh, 1:ntd] <= 0.
     p_fc_E[1:nh, 1:ntd] >= 0.
     p_tes_ch[1:nh, 1:ntd] <= 0.
     p_tes_dch[1:nh, 1:ntd] >= 0.
     p_heater_E[1:nh, 1:ntd] <= 0.
     p_g_out[1:nh, 1:ntd] <= 0.
     p_g_in[1:nh, 1:ntd] >= 0.
     # Investment decisions variables
     0 <= r_pv <= 1000
     0 <= r_liion <= 1000
     0 <= r_tank <= 50000
     0 <= r_elyz <= 50
     0 <= r_fc <= 50
     0 <= r_tes <= 1000
     # Operation state variables
     soc_liion[1:nh+1, 1:nd]
     soc_h2[1:nh+1, 1:nd]
     soc_tes[1:nh+1, 1:nd]
    end)

    # Constraints
    @constraints(m, begin
    # Power bounds
    [h in 1:nh, td in 1:ntd], p_liion_dch[h,td] <= liion.α_p_dch * r_liion
    [h in 1:nh, td in 1:ntd], p_liion_ch[h,td] >= -liion.α_p_ch * r_liion
    [h in 1:nh, td in 1:ntd], p_tes_dch[h,td] <= tes.α_p_dch * r_tes
    [h in 1:nh, td in 1:ntd], p_tes_ch[h,td] >= -tes.α_p_ch * r_tes
    [h in 1:nh, td in 1:ntd], p_fc_E[h,td] <= r_fc
    [h in 1:nh, td in 1:ntd], p_elyz_E[h,td] >= -r_elyz
    [h in 1:nh, td in 1:ntd], p_heater_E[h,td] >= -heater.powerMax[1]
    [h in 1:nh, td in 1:ntd], p_g_in[h,td] <= (1. - grid.τ_power) * maximum(ω_optim.clusters.ld_E .+ ω_optim.clusters.ld_H ./ heater.η_E_H) # seems to be faster in vectorized form
    # SoC bounds
    [h in 1:nh+1, d in 1:nd], soc_liion[h,d] <= liion.α_soc_max * r_liion
    [h in 1:nh+1, d in 1:nd], soc_liion[h,d] >= liion.α_soc_min * r_liion
    [h in 1:nh+1, d in 1:nd], soc_h2[h,d] <= h2tank.α_soc_max * r_tank
    [h in 1:nh+1, d in 1:nd], soc_h2[h,d] >= h2tank.α_soc_min * r_tank
    [h in 1:nh+1, d in 1:nd], soc_tes[h,d] <= tes.α_soc_max * r_tes
    [h in 1:nh+1, d in 1:nd], soc_tes[h,d] >= tes.α_soc_min * r_tes
    # State dynamics
    [h in 1:nh, d in 1:nd], soc_liion[h+1,d] == soc_liion[h,d] * (1. - liion.η_self * parameters.Δh) - (p_liion_ch[h,ω_td.σ[d]] * liion.η_ch + p_liion_dch[h,ω_td.σ[d]] / liion.η_dch) * parameters.Δh
    [h in 1:nh, d in 1:nd], soc_h2[h+1,d] == soc_h2[h,d] * (1. - h2tank.η_self * parameters.Δh) - (p_elyz_E[h,ω_td.σ[d]] * elyz.η_E_H2 + p_fc_E[h,ω_td.σ[d]] / fc.η_H2_E) * parameters.Δh
    [h in 1:nh, d in 1:nd], soc_tes[h+1,d] == soc_tes[h,d] * (1. - tes.η_self * parameters.Δh) - (p_tes_ch[h,ω_td.σ[d]] * tes.η_ch + p_tes_dch[h,ω_td.σ[d]] / tes.η_dch) * parameters.Δh
    # Sate continuity between days
    [d in 1:nd-1], soc_liion[1,d+1] == soc_liion[nh+1,d] * (1. - liion.η_self * parameters.Δh) - (p_liion_ch[1,ω_td.σ[d+1]] * liion.η_ch + p_liion_dch[1,ω_td.σ[d+1]] / liion.η_dch) * parameters.Δh
    [d in 1:nd-1], soc_h2[1,d+1] == soc_h2[nh+1,d] * (1. - h2tank.η_self * parameters.Δh) - (p_elyz_E[1,ω_td.σ[d+1]] * elyz.η_E_H2 + p_fc_E[1,ω_td.σ[d+1]] / fc.η_H2_E) * parameters.Δh
    [d in 1:nd-1], soc_tes[1,d+1] == soc_tes[nh+1,d] * (1. - tes.η_self * parameters.Δh) - (p_tes_ch[1,ω_td.σ[d+1]] * tes.η_ch + p_tes_dch[1,ω_td.σ[d+1]] / tes.η_dch) * parameters.Δh
    # Power balance
    [h in 1:nh, td in 1:ntd], ω_optim.clusters.ld_E[h,td] - r_pv * ω_optim.clusters.pv_E[h,td] <= p_g_out[h,td] + p_g_in[h,td] + p_liion_ch[h,td] + p_liion_dch[h,td] + p_elyz_E[h,td] + p_fc_E[h,td] + p_heater_E[h,td]
    [h in 1:nh, td in 1:ntd], ω_optim.clusters.ld_H[h,td] <= - elyz.η_E_H * p_elyz_E[h,td] +  fc.η_H2_H / fc.η_H2_E * p_fc_E[h,td] -  heater.η_E_H * p_heater_E[h,td] + p_tes_ch[h,td] + p_tes_dch[h,td]
    # Self-sufficiency constraint
    self_constraint, sum(ω_optim.nby[td] * sum(p_g_in[h,td] for h in 1:nh) for td in 1:ntd) <= (1. - grid.τ_energy) * sum(ω_optim.nby[td] * sum(ω_optim.clusters.ld_E[h,td] + ω_optim.clusters.ld_H[h,td] / heater.η_E_H for h in 1:nh) for td in 1:ntd)
    # Initial and final conditions
    soc_liion[1,1] == liion.soc[1,1] * r_liion
    soc_h2[1,1] == h2tank.soc[1,1] * r_tank
    soc_tes[1,1] == tes.soc[1,1] * r_tes
    soc_liion[nh,nd] >= soc_liion[1,1]
    soc_h2[nh,nd] >= soc_h2[1,1]
    soc_tes[nh,nd] >= soc_tes[1,1]
    end)

    # CAPEX
    capex = @expression(m, Γ_pv * ω_optim.clusters.C_pv[1] * r_pv +
    Γ_liion * ω_optim.clusters.C_liion[1] * r_liion +
    Γ_tank * ω_optim.clusters.C_tank[1] * r_tank +
    Γ_elyz * ω_optim.clusters.C_elyz[1] * r_elyz +
    Γ_fc * ω_optim.clusters.C_fc[1] * r_fc +
    Γ_tes * ω_optim.clusters.C_tes[1] * r_tes)

    # OPEX
    opex = @expression(m, sum(ω_optim.nby[td] * sum((p_g_in[h,td] * ω_optim.clusters.C_grid_in[h] + p_g_out[h,td] * ω_optim.clusters.C_grid_out[h]) * parameters.Δh  for h in 1:nh) for td in 1:ntd))

    # Objective
    @objective(m, Min, capex + opex)

    return m
end

#### Offline functions ####
# Multi-energy
function offline_optimization(ld::Load, pv::Source, liion::Liion, h2tank::H2Tank,
   elyz::Electrolyzer, fc::FuelCell, tes::ThermalSto, heater::Heater,
   designer::TypicalDayEACDesigner, grid::Grid, ω_optim::Scenarios,
   parameters::NamedTuple)

   # Parameters
   ny = size(ld.power_E,2) # number of simulation years
   ns = size(ld.power_E,3) # number of scenarios

   # Scenario reduction from the optimization scenario pool
   ω_eac = scenarios_reduction(ω_optim, "eac")

   # Clustering typical days
   ω_td = clustering_typical_days(ω_eac, designer.ntd)

   # Initialize model
   designer.model = eac_milp_model(ld, pv, liion, h2tank, elyz, fc, tes, heater, designer, grid, ω_td, parameters)

   # Compute both investment and operation decisions
   optimize!(designer.model)

    # Formatting variables to simulation
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
# Multi-energy
function compute_investment_decisions(y::Int64, s::Int64, ld::Load, pv::Source,
    liion::Liion, h2tank::H2Tank, elyz::Electrolyzer, fc::FuelCell, tes::ThermalSto,
    heater::Heater, designer::TypicalDayEACDesigner, ω_optim::Scenarios, parameters::NamedTuple)
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
