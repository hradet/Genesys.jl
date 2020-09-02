#=
    Designer based on the equivalent annual cost (EAC)
=#

mutable struct EACDesigner <: AbstractDesigner
    u::NamedTuple
    horizon::Int64
    model::JuMP.Model
    EACDesigner() = new()
end

#### Models ####
# Simple
function eac_milp_model(ld::Load, pv::Source, liion::Liion,
  designer::EACDesigner, grid::Grid, ω_optim::Scenarios,
  parameters::NamedTuple)

    # Parameters
    Γ_liion = (parameters.τ * (parameters.τ + 1.) ^ liion.lifetime) / ((parameters.τ + 1.) ^ liion.lifetime - 1.)
    Γ_pv = (parameters.τ * (parameters.τ + 1.) ^ pv.lifetime) / ((parameters.τ + 1.) ^ pv.lifetime - 1.)

    # Sets
    nh = size(ω_optim.ld_E,1) # Number of hours

    # Model definition
    m = Model(CPLEX.Optimizer)

    # Variables
    @variables(m, begin
    # Operation decision variables
     p_liion_ch[1:nh] <= 0.
     p_liion_dch[1:nh] >= 0.
     p_g_out[1:nh] <= 0.
     p_g_in[1:nh] >= 0.
     # Investment decision variables
     0 <= r_liion <= 1000
     0 <= r_pv <= 1000
     # Operation states variables
     soc_liion[1:nh+1]
    end)

    # Constraints
    @constraints(m, begin
    # Power bounds
    [h in 1:nh], p_liion_dch[h] <= liion.α_p_dch * r_liion
    [h in 1:nh], p_liion_ch[h] >= -liion.α_p_ch * r_liion
    [h in 1:nh], p_g_in[h] <= (1. - grid.τ_power) * maximum(ω_optim.ld_E)
    # SoC bounds
    [h in 1:nh+1], soc_liion[h] <= liion.α_soc_max * r_liion
    [h in 1:nh+1], soc_liion[h] >= liion.α_soc_min * r_liion
    # Investment bounds
    r_liion <= 1000.
    r_pv <= 1000.
    # State dynamic
    [h in 1:nh], soc_liion[h+1] == soc_liion[h] * (1 - liion.η_self * parameters.Δh) - (p_liion_ch[h] * liion.η_ch + p_liion_dch[h] / liion.η_dch) * parameters.Δh
    # Power balance
    [h in 1:nh], ω_optim.ld_E[h] - r_pv * ω_optim.pv_E[h] <= p_g_out[h] + p_g_in[h] + p_liion_ch[h] + p_liion_dch[h]
    # Self-sufficiency constraint
    self_constraint, sum(p_g_in[h] for h in 1:nh) <= (1. - grid.τ_energy) * sum(ω_optim.ld_E[h] for h in 1:nh)
    # Initial and final conditions
    soc_liion[1] == liion.soc[1] * r_liion
    soc_liion[end] == soc_liion[1]
    end)

    # CAPEX
    capex = @expression(m, Γ_pv * ω_optim.C_pv[1] * r_pv + Γ_liion * ω_optim.C_liion[1] * r_liion)

    # OPEX
    opex = @expression(m, sum((p_g_in[h] * ω_optim.C_grid_in[h] + p_g_out[h] * ω_optim.C_grid_out[h]) * parameters.Δh  for h in 1:nh))

    # Objective
    @objective(m, Min, capex + opex)

    return m
end
# Multi-energy
function eac_milp_model(ld::Load, pv::Source, liion::Liion, h2tank::H2Tank,
   elyz::Electrolyzer, fc::FuelCell, tes::ThermalSto, heater::Heater,
    designer::EACDesigner, grid::Grid, ω_optim::Scenarios,
    parameters::NamedTuple)
    # Parameters
    Γ_pv = (parameters.τ * (parameters.τ + 1.) ^ pv.lifetime) / ((parameters.τ + 1.) ^ pv.lifetime - 1.)
    Γ_liion = (parameters.τ * (parameters.τ + 1.) ^ liion.lifetime) / ((parameters.τ + 1.) ^ liion.lifetime - 1.)
    Γ_tank = (parameters.τ * (parameters.τ + 1.) ^ h2tank.lifetime) / ((parameters.τ + 1.) ^ h2tank.lifetime - 1.)
    Γ_elyz = (parameters.τ * (parameters.τ + 1.) ^ elyz.lifetime) / ((parameters.τ + 1.) ^ elyz.lifetime - 1.)
    Γ_fc = (parameters.τ * (parameters.τ + 1.) ^ fc.lifetime) / ((parameters.τ + 1.) ^ fc.lifetime - 1.)
    Γ_tes = (parameters.τ * (parameters.τ + 1.) ^ tes.lifetime) / ((parameters.τ + 1.) ^ tes.lifetime - 1.)

    # Sets
    nh = size(ω_optim.ld_E,1) # Number of hours

    # Model definition
    m = Model(CPLEX.Optimizer)

    # Variables
    @variables(m, begin
    # Operation decisions variables
    p_liion_ch[1:nh] <= 0.
    p_liion_dch[1:nh] >= 0.
    p_elyz_E[1:nh] <= 0.
    p_fc_E[1:nh] >= 0.
    p_tes_ch[1:nh] <= 0.
    p_tes_dch[1:nh] >= 0.
    p_heater_E[1:nh] <= 0.
    p_g_out[1:nh] <= 0.
    p_g_in[1:nh] >= 0.
    # Investment decisions variables
    0 <= r_pv <= 1000
    0 <= r_liion <= 1000
    0 <= r_tank <= 50000
    0 <= r_elyz <= 50
    0 <= r_fc <= 50
    0 <= r_tes <= 1000
    # Operation state variables
    soc_liion[1:nh+1]
    soc_h2[1:nh+1]
    soc_tes[1:nh+1]
    end)

    # Constraints
    @constraints(m, begin
    # Power bounds
    [h in 1:nh], p_liion_dch[h] <= liion.α_p_dch * r_liion
    [h in 1:nh], p_liion_ch[h] >= -liion.α_p_ch * r_liion
    [h in 1:nh], p_tes_dch[h] <= tes.α_p_dch * r_tes
    [h in 1:nh], p_tes_ch[h] >= -tes.α_p_ch * r_tes
    [h in 1:nh], p_fc_E[h] <= r_fc
    [h in 1:nh], p_elyz_E[h] >= -r_elyz
    [h in 1:nh], p_heater_E[h] >= -heater.powerMax[1]
    [h in 1:nh], p_g_in[h] <= (1. - grid.τ_power) * maximum(ω_optim.ld_E .+ ω_optim.ld_H ./ heater.η_E_H) # seems to be faster in vectorized form
    # SoC bounds
    [h in 1:nh+1], soc_liion[h] <= liion.α_soc_max * r_liion
    [h in 1:nh+1], soc_liion[h] >= liion.α_soc_min * r_liion
    [h in 1:nh+1], soc_h2[h] <= h2tank.α_soc_max * r_tank
    [h in 1:nh+1], soc_h2[h] >= h2tank.α_soc_min * r_tank
    [h in 1:nh+1], soc_tes[h] <= tes.α_soc_max * r_tes
    [h in 1:nh+1], soc_tes[h] >= tes.α_soc_min * r_tes
    # State dynamics
    [h in 1:nh], soc_liion[h+1] == soc_liion[h] * (1. - liion.η_self * parameters.Δh) - (p_liion_ch[h] * liion.η_ch + p_liion_dch[h] / liion.η_dch) * parameters.Δh
    [h in 1:nh], soc_h2[h+1] == soc_h2[h] * (1. - h2tank.η_self * parameters.Δh) - (p_elyz_E[h] * elyz.η_E_H2 + p_fc_E[h] / fc.η_H2_E) * parameters.Δh
    [h in 1:nh], soc_tes[h+1] == soc_tes[h] * (1. - tes.η_self * parameters.Δh) - (p_tes_ch[h] * tes.η_ch + p_tes_dch[h] / tes.η_dch) * parameters.Δh
    # Power balance
    [h in 1:nh], ω_optim.ld_E[h] - r_pv * ω_optim.pv_E[h] <= p_g_out[h] + p_g_in[h] + p_liion_ch[h] + p_liion_dch[h] + p_elyz_E[h] + p_fc_E[h] + p_heater_E[h]
    [h in 1:nh], ω_optim.ld_H[h] <= - elyz.η_E_H * p_elyz_E[h] +  fc.η_H2_H / fc.η_H2_E * p_fc_E[h] -  heater.η_E_H * p_heater_E[h] + p_tes_ch[h] + p_tes_dch[h]
    # Self-sufficiency constraint
    self_constraint, sum(p_g_in[h] for h in 1:nh) <= (1. - grid.τ_energy) * sum(ω_optim.ld_E[h] .+ ω_optim.ld_H[h] ./ heater.η_E_H for h in 1:nh)
    # Initial and final conditions
    soc_liion[1] == liion.soc[1,1] * r_liion
    soc_h2[1] == h2tank.soc[1,1] * r_tank
    soc_tes[1] == tes.soc[1,1] * r_tes
    soc_liion[nh] >= soc_liion[1]
    soc_h2[nh] >= soc_h2[1]
    soc_tes[nh] >= soc_tes[1]
    end)

    # CAPEX
    capex = @expression(m, Γ_pv * ω_optim.C_pv[1] * r_pv +
    Γ_liion * ω_optim.C_liion[1] * r_liion +
    Γ_tank * ω_optim.C_tank[1] * r_tank +
    Γ_elyz * ω_optim.C_elyz[1] * r_elyz +
    Γ_fc * ω_optim.C_fc[1] * r_fc +
    Γ_tes * ω_optim.C_tes[1] * r_tes)

    # OPEX
    opex = @expression(m, sum((p_g_in[h] * ω_optim.C_grid_in[h] + p_g_out[h] * ω_optim.C_grid_out[h]) * parameters.Δh  for h in 1:nh))

    # Objective
    @objective(m, Min, capex + opex)

    return m
end
# Multi-energy with clustering TODO : pb avec les scenarios !!!
function eac_milp_model_clustering(ld::Load, pv::Source, liion::Liion, h2tank::H2Tank,
   elyz::Electrolyzer, fc::FuelCell, tes::ThermalSto, heater::Heater,
    designer::EACDesigner, grid::Grid, ω_td::ClusteredScenarios,
    parameters::NamedTuple)
    # Parameters
    Γ_pv = (parameters.τ * (parameters.τ + 1.) ^ pv.lifetime) / ((parameters.τ + 1.) ^ pv.lifetime - 1.)
    Γ_liion = (parameters.τ * (parameters.τ + 1.) ^ liion.lifetime) / ((parameters.τ + 1.) ^ liion.lifetime - 1.)
    Γ_tank = (parameters.τ * (parameters.τ + 1.) ^ h2tank.lifetime) / ((parameters.τ + 1.) ^ h2tank.lifetime - 1.)
    Γ_elyz = (parameters.τ * (parameters.τ + 1.) ^ elyz.lifetime) / ((parameters.τ + 1.) ^ elyz.lifetime - 1.)
    Γ_fc = (parameters.τ * (parameters.τ + 1.) ^ fc.lifetime) / ((parameters.τ + 1.) ^ fc.lifetime - 1.)
    Γ_tes = (parameters.τ * (parameters.τ + 1.) ^ tes.lifetime) / ((parameters.τ + 1.) ^ tes.lifetime - 1.)

    # Sets
    nh = size(ω_td.clusters.ld_E,1) # number of hours by td (=24)
    ntd = size(ω_td.clusters.ld_E,2) # number of td
    nd = size(ω_td.σ,1) # number of days over the year (length of the sequence)

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
    [h in 1:nh, td in 1:ntd], p_g_in[h,td] <= (1. - grid.τ_power) * maximum(ω_optim.ld_E .+ ω_optim.ld_H ./ heater.η_E_H) # seems to be faster in vectorized form
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
    [h in 1:nh, td in 1:ntd], ω_td.clusters.ld_E[h,td] - r_pv * ω_optim.pv_E[h,td] <= p_g_out[h,td] + p_g_in[h,td] + p_liion_ch[h,td] + p_liion_dch[h,td] + p_elyz_E[h,td] + p_fc_E[h,td] + p_heater_E[h,td]
    [h in 1:nh, td in 1:ntd], ω_td.clusters.ld_H[h,td] <= - elyz.η_E_H * p_elyz_E[h,td] +  fc.η_H2_H / fc.η_H2_E * p_fc_E[h,td] -  heater.η_E_H * p_heater_E[h,td] + p_tes_ch[h,td] + p_tes_dch[h,td]
    # Self-sufficiency constraint
    self_constraint, sum(ω_td.nby[td] * sum(p_g_in[h,td] for h in 1:nh) for td in 1:ntd) <= (1. - grid.τ_energy) * sum(ω_td.nby[td] * sum(ω_td.clusters.ld_E[h,td] + ω_td.clusters.ld_H[h,td] / heater.η_E_H for h in 1:nh) for td in 1:ntd)
    # Initial and final conditions
    soc_liion[1,1] == liion.soc[1,1] * r_liion
    soc_h2[1,1] == h2tank.soc[1,1] * r_tank
    soc_tes[1,1] == tes.soc[1,1] * r_tes
    soc_liion[nh,nd] >= soc_liion[1,1]
    soc_h2[nh,nd] >= soc_h2[1,1]
    soc_tes[nh,nd] >= soc_tes[1,1]
    end)

    # CAPEX
    capex = @expression(m, Γ_pv * ω_optim.C_pv[1] * r_pv +
    Γ_liion * ω_optim.C_liion[1] * r_liion +
    Γ_tank * ω_optim.C_tank[1] * r_tank +
    Γ_elyz * ω_optim.C_elyz[1] * r_elyz +
    Γ_fc * ω_optim.C_fc[1] * r_fc +
    Γ_tes * ω_optim.C_tes[1] * r_tes)

    # OPEX
    opex = @expression(m, sum(ω_td.nby[td] * sum((p_g_in[h,td] * ω_td.clusters.C_grid_in[h] + p_g_out[h,td] * ω_td.clusters.C_grid_out[h]) * parameters.Δh  for h in 1:nh) for td in 1:ntd))

    # Objective
    @objective(m, Min, capex + opex)

    return m
end

#### Offline functions ####
# Simple
function offline_optimization(ld::Load, pv::Source, liion::Liion,
     designer::EACDesigner, grid::Grid, ω_optim::Scenarios, parameters::NamedTuple)
     # Parameters
     ny = size(ld.power_E,2) # number of simulation years
     ns = size(ld.power_E,3) # number of scenarios

     # Scenario reduction from the optimization scenario pool
     ω_milp = scenario_reduction(ω_optim, mode = "oneyear")

     # Initialize model
     designer.model = eac_milp_model(ld, pv, liion, designer, grid, ω_milp, parameters)

     # Compute investment decisions
     optimize!(designer.model)

     # Formatting variables to simulation
     designer.u = (
     u_liion = repeat(vcat(value.(designer.model[:r_liion]), zeros(ny-1,1)), 1, ns)
     u_pv = repeat(vcat(value.(designer.model[:r_pv]), zeros(ny-1,1)), 1, ns)
     )
end
# Multi-energy
function offline_optimization(ld::Load, pv::Source, liion::Liion, h2tank::H2Tank,
   elyz::Electrolyzer, fc::FuelCell, tes::ThermalSto, heater::Heater,
   designer::EACDesigner, grid::Grid, ω_optim::Scenarios,
   parameters::NamedTuple; mode="reference")

   # Parameters
   ny = size(ld.power_E,2) # number of simulation years
   ns = size(ld.power_E,3) # number of scenarios

  if mode == "reference"
    # Scenario reduction from the optimization scenario pool
    ω_milp = scenario_reduction(ω_optim, mode = "oneyear")

    # Initialize model
    designer.model = controller.model = onestage_milp_model(
    ld, pv, liion, h2tank, elyz, fc, tes, heater, designer, grid, ω_milp, parameters)

    # Compute investment decisions
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
   elseif mode == "clustering"
     # Selection of one-year scenarios from the optimization dataset
     ω_milp = scenario_reduction(ω_optim, mode = "typicaldays")

     # Initialize model
     designer.model = eac_milp_model_clustering(
     ld, pv, liion, h2tank, elyz, fc, tes, heater, designer, grid, ω_milp, parameters)

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
   else
     println("Unkwonw mode... Please chose between 'reference', or 'clustering' " )
   end
end

#### Online functions ####
# Simple
function compute_investment_decisions(y::Int64, s::Int64, ld::Load, pv::Source,
    liion::Liion, grid::Grid, designer::EACDesigner, ω_optim::Scenarios, parameters::NamedTuple)
    ϵ = 0.1
    if liion.soh[end,y,s] < ϵ
        designer.u.u_liion[y,s] = liion.Erated[y,s]
    end
end
# Multi-energy
function compute_investment_decisions(y::Int64, s::Int64, ld::Load, pv::Source,
    liion::Liion, h2tank::H2Tank, elyz::Electrolyzer, fc::FuelCell, tes::ThermalSto,
    heater::Heater, designer::EACDesigner, ω_optim::Scenarios, parameters::NamedTuple)
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
