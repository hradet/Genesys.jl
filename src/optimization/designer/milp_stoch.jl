#=
    Designer based on the equivalent annual cost (EAC) with multiple scenarios
=#

mutable struct MILPStochOptions
  solver
  risk
  scenario_reduction::String # scenario reduction technique "manual" or "auto"
  share_constraint::Bool
  reopt::Bool
  s::Int64 # scenario index for the manual reduction technique
  range_y::UnitRange{Int64} # year range for the manual reduction technique

  MILPStochOptions(; solver = CPLEX,
                    risk = "esperance",
                    scenario_reduction = "manual",
                    share_constraint = true,
                    reopt=false,
                    s = 1, range_y = 1:20) =
                    new(solver, risk, scenario_reduction, share_constraint, reopt, s, range_y)
end

mutable struct MILPStoch <: AbstractOneStageStochasticDesigner
    options::MILPStochOptions
    u::NamedTuple
    model::JuMP.Model
    history::AbstractScenarios

    MILPStoch(; options = MILPStochOptions()) = new(options)
end

### Models
function build_model(des::DistributedEnergySystem, designer::MILPStoch, ω_optim::AbstractScenarios)

    # Sets
    nh = size(ω_optim.values.ld_E,1) # Number of hours
    ns = size(ω_optim.values.ld_E,2) # Number of scenarios

    # Model definition
    m = Model(designer.options.solver.Optimizer)

    # Initialize expressions
    isa(des.ld_E, Load) ? power_balance_E = AffExpr.(- ω_optim.values.ld_E) : nothing
    isa(des.ld_H, Load) ? power_balance_H = AffExpr.(- ω_optim.values.ld_H) : nothing
    power_balance_H2 = AffExpr.(zeros(nh, ns))
    capex = AffExpr(0.)

    # Build model
    if isa(des.liion, Liion)
        # Parameters
        Γ_liion = (des.parameters.τ * (des.parameters.τ + 1.) ^ des.liion.lifetime) / ((des.parameters.τ + 1.) ^ des.liion.lifetime - 1.)
        # Variables
        @variables(m, begin
        # Operation decisions variables
        p_liion_ch[1:nh, 1:ns] <= 0.
        p_liion_dch[1:nh, 1:ns] >= 0.
        # Investment decisions variables
        0 <= r_liion <= 1000
        # Operation state variables
        soc_liion[1:nh+1, 1:ns]
        end)
        # Constraints
        @constraints(m, begin
        # Power bounds
        [h in 1:nh, s in 1:ns], p_liion_dch[h,s] <= des.liion.α_p_dch * r_liion
        [h in 1:nh, s in 1:ns], p_liion_ch[h,s] >= -des.liion.α_p_ch * r_liion
        # SoC bounds
        [h in 1:nh+1, s in 1:ns], soc_liion[h,s] <= des.liion.α_soc_max * r_liion
        [h in 1:nh+1, s in 1:ns], soc_liion[h,s] >= des.liion.α_soc_min * r_liion
        # State dynamics
        [h in 1:nh, s in 1:ns], soc_liion[h+1,s] == soc_liion[h,s] * (1. - des.liion.η_self * des.parameters.Δh) - (p_liion_ch[h,s] * des.liion.η_ch + p_liion_dch[h,s] / des.liion.η_dch) * des.parameters.Δh
        # Initial and final conditions
        [s in 1:ns], soc_liion[1,s] == des.liion.soc_ini * r_liion
        [s in 1:ns], soc_liion[nh,s] >= soc_liion[1,s]
        end)
        # Power balance
        isa(des.ld_E, Load) ? add_to_expression!.(power_balance_E, p_liion_ch .+ p_liion_dch) : nothing
        # CAPEX
        add_to_expression!(capex, Γ_liion * ω_optim.values.C_liion[1] * r_liion)
    end

    if isa(des.tes, ThermalSto)
        # Parameters
        Γ_tes = (des.parameters.τ * (des.parameters.τ + 1.) ^ des.tes.lifetime) / ((des.parameters.τ + 1.) ^ des.tes.lifetime - 1.)
        # Variables
        @variables(m, begin
        # Operation decisions variables
        p_tes_ch[1:nh, 1:ns] <= 0.
        p_tes_dch[1:nh, 1:ns] >= 0.
        # Investment decisions variables
        0 <= r_tes <= 1000
        # Operation state variables
        soc_tes[1:nh+1, 1:ns]
        end)
        # Constraints
        @constraints(m, begin
        # Power bounds
        [h in 1:nh, s in 1:ns], p_tes_dch[h,s] <= des.tes.α_p_dch * r_tes
        [h in 1:nh, s in 1:ns], p_tes_ch[h,s] >= -des.tes.α_p_ch * r_tes
        # SoC bounds
        [h in 1:nh+1, s in 1:ns], soc_tes[h,s] <= des.tes.α_soc_max * r_tes
        [h in 1:nh+1, s in 1:ns], soc_tes[h,s] >= des.tes.α_soc_min * r_tes
        # State dynamics
        [h in 1:nh, s in 1:ns], soc_tes[h+1,s] == soc_tes[h,s] * (1. - des.tes.η_self * des.parameters.Δh) - (p_tes_ch[h,s] * des.tes.η_ch + p_tes_dch[h,s] / des.tes.η_dch) * des.parameters.Δh
        # Initial and final conditions
        [s in 1:ns], soc_tes[1,s] == des.tes.soc_ini * r_tes
        [s in 1:ns], soc_tes[nh,s] >= soc_tes[1,s]
        end)
        # Power balance
        isa(des.ld_H, Load) ? add_to_expression!.(power_balance_H, p_tes_ch .+ p_tes_dch) : nothing
        # CAPEX
        add_to_expression!(capex, Γ_tes * ω_optim.values.C_tes[1] * r_tes)
    end

    if isa(des.h2tank, H2Tank)
        # Parameters
        Γ_h2tank = (des.parameters.τ * (des.parameters.τ + 1.) ^ des.h2tank.lifetime) / ((des.parameters.τ + 1.) ^ des.h2tank.lifetime - 1.)
        # Variables
        @variables(m, begin
        # Operation decisions variables
        p_h2tank_ch[1:nh, 1:ns] <= 0.
        p_h2tank_dch[1:nh, 1:ns] >= 0.
        # Investment decisions variables
        0 <= r_h2tank <= 50000
        # Operation state variables
        soc_h2tank[1:nh+1, 1:ns]
        end)
        # Constraints
        @constraints(m, begin
        # Power bounds
        [h in 1:nh, s in 1:ns], p_h2tank_dch[h,s] <= des.h2tank.α_p_dch * r_h2tank
        [h in 1:nh, s in 1:ns], p_h2tank_ch[h,s] >= -des.h2tank.α_p_ch * r_h2tank
        # SoC bounds
        [h in 1:nh+1, s in 1:ns], soc_h2tank[h,s] <= des.h2tank.α_soc_max * r_h2tank
        [h in 1:nh+1, s in 1:ns], soc_h2tank[h,s] >= des.h2tank.α_soc_min * r_h2tank
        # State dynamics
        [h in 1:nh, s in 1:ns], soc_h2tank[h+1,s] == soc_h2tank[h,s] * (1. - des.h2tank.η_self * des.parameters.Δh) - (p_h2tank_ch[h,s] * des.h2tank.η_ch + p_h2tank_dch[h,s] / des.h2tank.η_dch) * des.parameters.Δh
        # Initial and final conditions
        [s in 1:ns], soc_h2tank[1,s] == des.h2tank.soc_ini * r_h2tank
        [s in 1:ns], soc_h2tank[nh,s] >= soc_h2tank[1,s]
        end)
        # Power balances
        add_to_expression!.(power_balance_H2, p_h2tank_ch .+ p_h2tank_dch)
        # CAPEX
        add_to_expression!.(capex, Γ_h2tank * ω_optim.values.C_tank[1] * r_h2tank)
    end

    if isa(des.elyz, Electrolyzer)
        # Parameters
        Γ_elyz = (des.parameters.τ * (des.parameters.τ + 1.) ^ des.elyz.lifetime) / ((des.parameters.τ + 1.) ^ des.elyz.lifetime - 1.)
        # Variables
        @variables(m, begin
        # Operation decisions variables
        p_elyz_E[1:nh, 1:ns] <= 0.
        # Investment decisions variables
        0 <= r_elyz <= 50
        end)
        # Constraints
        @constraints(m, begin
        # Power bounds
        [h in 1:nh, s in 1:ns], p_elyz_E[h,s] >= -r_elyz
        end)
        # Power balance
        isa(des.ld_E, Load) ? add_to_expression!.(power_balance_E, p_elyz_E) : nothing
        isa(des.ld_H, Load) ? add_to_expression!.(power_balance_H, - des.elyz.η_E_H .* p_elyz_E) : nothing
        add_to_expression!.(power_balance_H2, - des.elyz.η_E_H2 * p_elyz_E)
        # CAPEX
        add_to_expression!(capex, Γ_elyz * ω_optim.values.C_elyz[1] * r_elyz)
    end

    if isa(des.fc, FuelCell)
        # Parameters
        Γ_fc = (des.parameters.τ * (des.parameters.τ + 1.) ^ des.fc.lifetime) / ((des.parameters.τ + 1.) ^ des.fc.lifetime - 1.)
        # Variables
        @variables(m, begin
        # Operation decisions variables
        p_fc_E[1:nh, 1:ns] >= 0.
        # Investment decisions variables
        0 <= r_fc <= 50
        end)
        # Constraints
        @constraints(m, begin
        # Power bounds
        [h in 1:nh, s in 1:ns], p_fc_E[h,s] <= r_fc
        end)
        # Power balance
        isa(des.ld_E, Load) ? add_to_expression!.(power_balance_E, p_fc_E) : nothing
        isa(des.ld_H, Load) ? add_to_expression!.(power_balance_H, des.fc.η_H2_H / des.fc.η_H2_E .* p_fc_E) : nothing
        add_to_expression!.(power_balance_H2, - p_fc_E / des.fc.η_H2_E)
        # CAPEX
        add_to_expression!(capex, Γ_fc * ω_optim.values.C_fc[1] * r_fc)
    end

    if isa(des.heater, Heater)
        # Variables
        @variables(m, begin
        # Operation decisions variables
        p_heater_E[1:nh, 1:ns] <= 0.
        end)
        # Constraints
        @constraints(m, begin
        # Power bounds
        [h in 1:nh, s in 1:ns], p_heater_E[h,s] >= -des.heater.powerMax[1]
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
        isa(des.ld_E, Load) ? add_to_expression!.(power_balance_E, r_pv .* ω_optim.values.pv_E) : nothing
        # CAPEX
        add_to_expression!(capex, Γ_pv * ω_optim.values.C_pv[1] * r_pv)
    end

    if isa(des.grid, Grid)
        # Variables
        @variables(m, begin
        # Operation decisions variables
        p_g_out[1:nh, 1:ns] <= 0.
        p_g_in[1:nh, 1:ns] >= 0.
        end)
        # Constraints
        @constraints(m, begin
        # Power bounds
        [h in 1:nh, s in 1:ns], p_g_in[h,s] <= des.grid.powerMax
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
        ld_tot = zeros(ns)
        isa(des.ld_E, Load) ? ld_tot .+= sum(ω_optim.values.ld_E[h,:] for h in 1:nh) : nothing
        isa(des.ld_H, Load) ? ld_tot .+= sum(ω_optim.values.ld_H[h,:] ./ des.heater.η_E_H for h in 1:nh) : nothing
        @constraint(m, self_constraint[s in 1:ns], sum(p_g_in[h,s] for h in 1:nh) <= (1. - des.parameters.τ_share) * ld_tot[s])
    end

    # OPEX
    if designer.options.risk == "esperance"
        opex = @expression(m, sum((p_g_in[h,s] * ω_optim.values.C_grid_in[h,s] + p_g_out[h,s] * ω_optim.values.C_grid_out[h,s]) * des.parameters.Δh  for h in 1:nh, s in 1:ns) / ns)
    elseif designer.options.risk == "cvar"
    else
      println("Unknown risk measure... Chose between 'esperance' or 'cvar' ")
    end

    # Objective
    @objective(m, Min, capex + opex)

    return m
end

### Offline
function initialize_designer!(des::DistributedEnergySystem, designer::MILPStoch, ω_optim::AbstractScenarios)

     # Save history for online optimization
     designer.history = ω_optim

     # Preallocate
     preallocate!(designer, des.parameters.ny, des.parameters.ns)

     return designer
end

### Online
function compute_investment_decisions!(y::Int64, s::Int64, des::DistributedEnergySystem, designer::MILPStoch)
    ϵ = 0.1

    if s == 1 && y == 1
        # Scenario reduction from the optimization scenario pool
        ω_eac_stoch = scenarios_reduction(designer, designer.history)

        # Initialize model
        designer.model = build_model(des, designer, ω_eac_stoch)

        # Compute investment decisions
        optimize!(designer.model)

        # Assign values
        isa(des.pv, Source) ? designer.u.pv[1,:] .= value(designer.model[:r_pv]) : nothing
        isa(des.liion, Liion) ? designer.u.liion[1,:] .= value(designer.model[:r_liion]) : nothing
        isa(des.h2tank, H2Tank) ? designer.u.h2tank[1,:] .= value(designer.model[:r_h2tank]) : nothing
        isa(des.elyz, Electrolyzer) ? designer.u.elyz[1,:] .= value(designer.model[:r_elyz]) : nothing
        isa(des.fc, FuelCell) ? designer.u.fc[1,:] .= value(designer.model[:r_fc]) : nothing
        isa(des.tes, ThermalSto) ? designer.u.tes[1,:] .= value(designer.model[:r_tes]) : nothing

    elseif designer.options.reopt
        # Do we need to reoptimize ?
        (isa(des.liion, Liion) && des.liion.soh[end,y,s] < ϵ) || (isa(des.elyz, Electrolyzer) && des.elyz.soh[end,y,s] < ϵ) || (isa(des.fc, FuelCell) && des.fc.soh[end,y,s] < ϵ) ? nothing : return

        # Update designer options with current year for scenario reduction
        designer.options.range_y = y:des.parameters.ny

        # Scenario reduction from the optimization scenario pool
        ω_eac_stoch = scenarios_reduction(designer, designer.history)

        # Build model
        designer.model = build_model(des, designer, ω_eac_stoch)

        # Fix variables
        isa(des.pv, Source) ? fix(designer.model[:r_pv], des.pv.powerMax[y,s], force=true) : nothing
        isa(des.h2tank, H2Tank) ? fix(designer.model[:r_h2tank], des.h2tank.Erated[y,s], force=true) : nothing
        isa(des.tes, ThermalSto) ? fix(designer.model[:r_tes], des.tes.Erated[y,s], force=true) : nothing
        isa(des.liion, Liion) && des.liion.soh[end,y,s] >= ϵ ? fix(designer.model[:r_liion], des.liion.Erated[y,s], force=true) : nothing
        isa(des.elyz, Electrolyzer) && des.elyz.soh[end,y,s] >= ϵ ? fix(designer.model[:r_elyz], des.elyz.powerMax[y,s], force=true) : nothing
        isa(des.fc, FuelCell) && des.fc.soh[end,y,s] >= ϵ ? fix(designer.model[:r_fc], des.fc.powerMax[y,s], force=true) : nothing

        # Compute investment decisions
        optimize!(designer.model)

        # Preallocate and assigned values
        isa(des.liion, Liion) && des.liion.soh[end,y,s] < ϵ ? designer.u.liion[y,s] = value(designer.model[:r_liion]) : nothing
        isa(des.elyz, Electrolyzer) && des.elyz.soh[end,y,s] < ϵ ? designer.u.elyz[y,s] = value(designer.model[:r_elyz]) : nothing
        isa(des.fc, FuelCell) && des.fc.soh[end,y,s] < ϵ ? designer.u.fc[y,s] = value(designer.model[:r_fc]) : nothing

    else
        isa(des.liion, Liion) && des.liion.soh[end,y,s] < ϵ ? designer.u.liion[y,s] = des.liion.Erated[y,s] : nothing
        isa(des.elyz, Electrolyzer) && des.elyz.soh[end,y,s] < ϵ ? designer.u.elyz[y,s] = des.elyz.powerMax[y,s] : nothing
        isa(des.fc, FuelCell) && des.fc.soh[end,y,s] < ϵ ? designer.u.fc[y,s] = des.fc.powerMax[y,s] : nothing
    end
end
