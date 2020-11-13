#=
    Equivalent annual cost (EAC) designer based on a single year
=#

mutable struct MILPOptions
  solver
  scenario_reduction::String # scenario reduction technique "manual" or "auto"
  share_constraint::Bool
  reopt::Bool
  s::Int64 # scenario index for the manual reduction technique
  y::Int64 # year index for the manual reduction technique

  MILPOptions(; solver = CPLEX,
               scenario_reduction = "manual",
               share_constraint = true,
               reopt = false,
               s = 1, y = 1) =
               new(solver, scenario_reduction, share_constraint, reopt, s, y)
end

mutable struct MILP <: AbstractDesigner
    options::MILPOptions
    u::NamedTuple
    model::JuMP.Model
    history::AbstractScenarios

    MILP(; options = MILPOptions()) = new(options)
end

### Model
function build_model(des::DistributedEnergySystem, designer::MILP, ω::AbstractScenarios)

    #TODO diviser en fonctions !! add_constraints!(model::JuMP.model, liion::Liion)
    # ou ecrire un seul probleme et fixer les variables de design à 0...

    # Sets
    nh = des.parameters.nh # Number of hours

    # Model definition
    m = Model(designer.options.solver.Optimizer)

    # Initialize expressions
    isa(des.ld_E, Load) ? power_balance_E = AffExpr.(- ω.ld_E.power) : nothing
    isa(des.ld_H, Load) ? power_balance_H = AffExpr.(- ω.ld_H.power) : nothing
    power_balance_H2 = AffExpr.(zeros(nh))
    capex = AffExpr(0.)

    # Build model
    if isa(des.liion, Liion)
        # Parameters
        Γ_liion = (des.parameters.τ * (des.parameters.τ + 1.) ^ des.liion.lifetime) / ((des.parameters.τ + 1.) ^ des.liion.lifetime - 1.)
        # Variables
        @variables(m, begin
        # Operation decisions variables
        p_liion_ch[1:nh] <= 0.
        p_liion_dch[1:nh] >= 0.
        # Investment decisions variables
        0 <= r_liion <= 1000
        # Operation state variables
        soc_liion[1:nh+1]
        end)
        # Constraints
        @constraints(m, begin
        # Power bounds
        [h in 1:nh], p_liion_dch[h] <= des.liion.α_p_dch * r_liion
        [h in 1:nh], p_liion_ch[h] >= -des.liion.α_p_ch * r_liion
        # SoC bounds
        [h in 1:nh+1], soc_liion[h] <= des.liion.α_soc_max * r_liion
        [h in 1:nh+1], soc_liion[h] >= des.liion.α_soc_min * r_liion
        # State dynamics
        [h in 1:nh], soc_liion[h+1] == soc_liion[h] * (1. - des.liion.η_self * des.parameters.Δh) - (p_liion_ch[h] * des.liion.η_ch + p_liion_dch[h] / des.liion.η_dch) * des.parameters.Δh
        # Initial and final conditions
        soc_liion[1] == des.liion.soc_ini * r_liion
        soc_liion[nh] >= soc_liion[1]
        end)
        # Power balance
        isa(des.ld_E, Load) ? add_to_expression!.(power_balance_E, p_liion_ch .+ p_liion_dch) : nothing
        # CAPEX
        add_to_expression!(capex, Γ_liion * ω.liion.cost[1] * r_liion)
    end

    if isa(des.tes, ThermalSto)
        # Parameters
        Γ_tes = (des.parameters.τ * (des.parameters.τ + 1.) ^ des.tes.lifetime) / ((des.parameters.τ + 1.) ^ des.tes.lifetime - 1.)
        # Variables
        @variables(m, begin
        # Operation decisions variables
        p_tes_ch[1:nh] <= 0.
        p_tes_dch[1:nh] >= 0.
        # Investment decisions variables
        0 <= r_tes <= 1000
        # Operation state variables
        soc_tes[1:nh+1]
        end)
        # Constraints
        @constraints(m, begin
        # Power bounds
        [h in 1:nh], p_tes_dch[h] <= des.tes.α_p_dch * r_tes
        [h in 1:nh], p_tes_ch[h] >= -des.tes.α_p_ch * r_tes
        # SoC bounds
        [h in 1:nh+1], soc_tes[h] <= des.tes.α_soc_max * r_tes
        [h in 1:nh+1], soc_tes[h] >= des.tes.α_soc_min * r_tes
        # State dynamics
        [h in 1:nh], soc_tes[h+1] == soc_tes[h] * (1. - des.tes.η_self * des.parameters.Δh) - (p_tes_ch[h] * des.tes.η_ch + p_tes_dch[h] / des.tes.η_dch) * des.parameters.Δh
        # Initial and final conditions
        soc_tes[1] == des.tes.soc_ini * r_tes
        soc_tes[nh] >= soc_tes[1]
        end)
        # Power balance
        isa(des.ld_H, Load) ? add_to_expression!.(power_balance_H, p_tes_ch .+ p_tes_dch) : nothing
        # CAPEX
        add_to_expression!(capex, Γ_tes * ω.tes.cost[1] * r_tes)
    end

    if isa(des.h2tank, H2Tank)
        # Parameters
        Γ_h2tank = (des.parameters.τ * (des.parameters.τ + 1.) ^ des.h2tank.lifetime) / ((des.parameters.τ + 1.) ^ des.h2tank.lifetime - 1.)
        # Variables
        @variables(m, begin
        # Operation decisions variables
        p_h2tank_ch[1:nh] <= 0.
        p_h2tank_dch[1:nh] >= 0.
        # Investment decisions variables
        0 <= r_h2tank <= 50000
        # Operation state variables
        soc_h2tank[1:nh+1]
        end)
        # Constraints
        @constraints(m, begin
        # Power bounds
        [h in 1:nh], p_h2tank_dch[h] <= des.h2tank.α_p_dch * r_h2tank
        [h in 1:nh], p_h2tank_ch[h] >= -des.h2tank.α_p_ch * r_h2tank
        # SoC bounds
        [h in 1:nh+1], soc_h2tank[h] <= des.h2tank.α_soc_max * r_h2tank
        [h in 1:nh+1], soc_h2tank[h] >= des.h2tank.α_soc_min * r_h2tank
        # State dynamics
        [h in 1:nh], soc_h2tank[h+1] == soc_h2tank[h] * (1. - des.h2tank.η_self * des.parameters.Δh) - (p_h2tank_ch[h] * des.h2tank.η_ch + p_h2tank_dch[h] / des.h2tank.η_dch) * des.parameters.Δh
        # Initial and final conditions
        soc_h2tank[1] == des.h2tank.soc_ini * r_h2tank
        soc_h2tank[nh] >= soc_h2tank[1]
        end)
        # Power balances
        add_to_expression!.(power_balance_H2, p_h2tank_ch .+ p_h2tank_dch)
        # CAPEX
        add_to_expression!.(capex, Γ_h2tank * ω.h2tank.cost[1] * r_h2tank)
    end

    if isa(des.elyz, Electrolyzer)
        # Parameters
        Γ_elyz = (des.parameters.τ * (des.parameters.τ + 1.) ^ des.elyz.lifetime) / ((des.parameters.τ + 1.) ^ des.elyz.lifetime - 1.)
        # Variables
        @variables(m, begin
        # Operation decisions variables
        p_elyz_E[1:nh] <= 0.
        # Investment decisions variables
        0 <= r_elyz <= 50
        end)
        # Constraints
        @constraints(m, begin
        # Power bounds
        [h in 1:nh], p_elyz_E[h] >= -r_elyz
        end)
        # Power balance
        isa(des.ld_E, Load) ? add_to_expression!.(power_balance_E, p_elyz_E) : nothing
        isa(des.ld_H, Load) ? add_to_expression!.(power_balance_H, - des.elyz.η_E_H .* p_elyz_E) : nothing
        add_to_expression!.(power_balance_H2, - des.elyz.η_E_H2 * p_elyz_E)
        # CAPEX
        add_to_expression!(capex, Γ_elyz * ω.elyz.cost[1] * r_elyz)
    end

    if isa(des.fc, FuelCell)
        # Parameters
        Γ_fc = (des.parameters.τ * (des.parameters.τ + 1.) ^ des.fc.lifetime) / ((des.parameters.τ + 1.) ^ des.fc.lifetime - 1.)
        # Variables
        @variables(m, begin
        # Operation decisions variables
        p_fc_E[1:nh] >= 0.
        # Investment decisions variables
        0 <= r_fc <= 50
        end)
        # Constraints
        @constraints(m, begin
        # Power bounds
        [h in 1:nh], p_fc_E[h] <= r_fc
        end)
        # Power balance
        isa(des.ld_E, Load) ? add_to_expression!.(power_balance_E, p_fc_E) : nothing
        isa(des.ld_H, Load) ? add_to_expression!.(power_balance_H, des.fc.η_H2_H / des.fc.η_H2_E .* p_fc_E) : nothing
        add_to_expression!.(power_balance_H2, - p_fc_E / des.fc.η_H2_E)
        # CAPEX
        add_to_expression!(capex, Γ_fc * ω.fc.cost[1] * r_fc)
    end

    if isa(des.heater, Heater)
        # Variables
        @variables(m, begin
        # Operation decisions variables
        p_heater_E[1:nh] <= 0.
        end)
        # Constraints
        @constraints(m, begin
        # Power bounds
        [h in 1:nh], p_heater_E[h] >= -des.heater.powerMax[1]
        end)
        # Power balance
        isa(des.ld_E, Load) ? add_to_expression!.(power_balance_E, p_heater_E) : nothing
        isa(des.ld_H, Load) ? add_to_expression!.(power_balance_H, - des.heater.η_E_H .* p_heater_E) : nothing
    end

    if isa(des.pv, Source)
        # Parameters
        Γ_pv = (des.parameters.τ * (des.parameters.τ + 1.) ^ des.pv.lifetime) / ((des.parameters.τ + 1.) ^ des.pv.lifetime - 1.)
        # Variables
        @variables(m, begin
        # Investment decisions variables
        0 <= r_pv <= 1000
        end)
        # Power balance
        isa(des.ld_E, Load) ? add_to_expression!.(power_balance_E, r_pv * ω.pv.power) : nothing
        # CAPEX
        add_to_expression!(capex, Γ_pv * ω.pv.cost[1] * r_pv)
    end

    if isa(des.grid, Grid)
        # Variables
        @variables(m, begin
        # Operation decisions variables
        p_g_out[1:nh] <= 0.
        p_g_in[1:nh] >= 0.
        end)
        # Constraints
        @constraints(m, begin
        # Power bounds
        [h in 1:nh], p_g_in[h] <= des.grid.powerMax
        end)
        # Power balance
        isa(des.ld_E, Load) ? add_to_expression!.(power_balance_E, p_g_in .+ p_g_out) : nothing
    end

    # Power balances
    isa(des.ld_E, Load) ? @constraint(m, power_balance_E .>= 0.) : nothing
    isa(des.ld_H, Load) ? @constraint(m, power_balance_H .>= 0.) : nothing
    @constraint(m, power_balance_H2 .== 0.)

    # Share of renewables constraint
    if designer.options.share_constraint
        sum_ld = (isa(des.ld_E, Load) ? sum(ω.ld_E.power[h] for h in 1:nh) : 0) + (isa(des.ld_H, Load) ? sum(ω.ld_H.power[h] ./ des.heater.η_E_H for h in 1:nh) : 0)
        @constraint(m, self_constraint, sum(p_g_in[h] for h in 1:nh) <= (1. - des.parameters.τ_share) * sum_ld)
    end

    # OPEX
    opex = @expression(m, sum((p_g_in[h] * ω.grid.cost_in[h] + p_g_out[h] * ω.grid.cost_out[h]) * des.parameters.Δh  for h in 1:nh))

    # Objective
    @objective(m, Min, capex + opex)

    return m
end

### Offline
function initialize_designer!(des::DistributedEnergySystem, designer::MILP, ω::AbstractScenarios)

   # Save history for online optimization
   designer.history = ω

   # Preallocate
   preallocate!(designer, des.parameters.ny, des.parameters.ns)

   return designer
end

### Online
function compute_investment_decisions!(y::Int64, s::Int64, des::DistributedEnergySystem, designer::MILP)
    # Parameters
    ϵ = 0.1

    if s == 1 && y == 1
        # Scenario reduction from the optimization scenario pool
        ω = scenarios_reduction(designer.history, 1:des.parameters.nh, 1, 1)

        # Build model
        designer.model = build_model(des, designer, ω)

        # Compute investment decisions
        optimize!(designer.model)

        # Assign values
        isa(des.pv, Source) ? designer.u.pv[1,:] .= value(designer.model[:r_pv]) : nothing
        isa(des.liion, Liion) ? designer.u.liion[1,:] .= value(designer.model[:r_liion]) : nothing
        isa(des.h2tank, H2Tank) ? designer.u.h2tank[1,:] .= value(designer.model[:r_h2tank]) : nothing
        isa(des.elyz, Electrolyzer) ? designer.u.elyz[1,:] .= value(designer.model[:r_elyz]) : nothing
        isa(des.fc, FuelCell) ? designer.u.fc[1,:] .= value(designer.model[:r_fc]) : nothing
        isa(des.tes, ThermalSto) ? designer.u.tes[1,:] .= value(designer.model[:r_tes]) : nothing

    elseif designer.options.reopt # reoptimization
        # Do we need to reoptimize ?
        (isa(des.liion, Liion) && des.liion.soh[end,y,s] < ϵ) || (isa(des.elyz, Electrolyzer) && des.elyz.soh[end,y,s] < ϵ) || (isa(des.fc, FuelCell) && des.fc.soh[end,y,s] < ϵ) ? nothing : return

        # Scenario reduction from the optimization scenario pool
        ω = scenarios_reduction(designer.history, 1:des.parameters.nh, y, 1)

        # Build model
        designer.model = build_model(des, designer, ω)

        # Fix variables
        isa(des.pv, Source) ? fix(designer.model[:r_pv], des.pv.powerMax[y,s], force=true) : nothing
        isa(des.h2tank, H2Tank) ? fix(designer.model[:r_h2tank], des.h2tank.Erated[y,s], force=true) : nothing
        isa(des.tes, ThermalSto) ? fix(designer.model[:r_tes], des.tes.Erated[y,s], force=true) : nothing
        isa(des.liion, Liion) && des.liion.soh[end,y,s] >= ϵ ? fix(designer.model[:r_liion], des.liion.Erated[y,s], force=true) : nothing
        isa(des.elyz, Electrolyzer) && des.elyz.soh[end,y,s] >= ϵ ? fix(designer.model[:r_elyz], des.elyz.powerMax[y,s], force=true) : nothing
        isa(des.fc, FuelCell) && des.fc.soh[end,y,s] >= ϵ ? fix(designer.model[:r_fc], des.fc.powerMax[y,s], force=true) : nothing

        # Compute investment decisions
        optimize!(designer.model)

        # Assign values
        isa(des.liion, Liion) && des.liion.soh[end,y,s] < ϵ ? designer.u.liion[y,s] = value(designer.model[:r_liion]) : nothing
        isa(des.elyz, Electrolyzer) && des.elyz.soh[end,y,s] < ϵ ? designer.u.elyz[y,s] = value(designer.model[:r_elyz]) : nothing
        isa(des.fc, FuelCell) && des.fc.soh[end,y,s] < ϵ ? designer.u.fc[y,s] = value(designer.model[:r_fc]) : nothing

    else
        isa(des.liion, Liion) && des.liion.soh[end,y,s] < ϵ ? designer.u.liion[y,s] = designer.u.liion[1,s] : nothing
        isa(des.elyz, Electrolyzer) && des.elyz.soh[end,y,s] < ϵ ? designer.u.elyz[y,s] = designer.u.elyz[1,s] : nothing
        isa(des.fc, FuelCell) && des.fc.soh[end,y,s] < ϵ ? designer.u.fc[y,s] = designer.u.fc[1,s] : nothing
    end
end
