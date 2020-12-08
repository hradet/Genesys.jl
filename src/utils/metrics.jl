#=
    This file includes all the funtions needed to compute the techno-economic
    indicators
 =#


 mutable struct Costs
      capex::Union{Array{Float64,1}, Array{Float64,2}}
      opex::Union{Array{Float64,1}, Array{Float64,2}}
      salvage::Union{Array{Float64,1}, Array{Float64,2}}
      cf::Union{Array{Float64,1}, Array{Float64,2}}
      cumulative_npv::Union{Array{Float64,1}, Array{Float64,2}}
      npv::Union{Float64, Array{Float64,1}}
 end

 # Compute costs
function Costs(des::DistributedEnergySystem, designer::AbstractDesigner)

     # Discount factor
     γ = 1. ./ (1. + des.parameters.τ) .^ range(0, length = des.parameters.ny, step = des.parameters.Δy)

    # Discounted capex
    capex = γ .* compute_capex(des, designer)

    # Discounted opex
    opex = γ .* (compute_baseline_cost(des) .- compute_grid_cost(des))

    # Discounted salvage value
    salvage = γ .* compute_salvage(des)

    # Discounted cash flow
    cf = - capex .+ opex .+ salvage

    # Discounted NPV each year
    cumulative_npv = cumsum(cf, dims = 1)

    # Discounted NPV
    npv = dropdims(sum(cf, dims=1), dims=1)

    return Costs(capex, opex, salvage, cf, cumulative_npv, npv)
end
# Compute costs for a given scenario s
function Costs(s::Int64, des::DistributedEnergySystem, designer::AbstractDesigner)

     # Discount factor
     γ = 1. ./ (1. + des.parameters.τ) .^ range(0, length = des.parameters.ny, step = des.parameters.Δy)

    # Discounted capex
    capex = γ .* compute_capex(s, des, designer)

    # Discounted opex
    opex = γ .* (compute_baseline_cost(s, des) .- compute_grid_cost(s, des))

    # Discounted salvage value
    salvage = γ .* compute_salvage(s, des)

    # Discounted cash flow
    cf = - capex .+ opex .+ salvage

    # Discounted NPV each year
    cumulative_npv = cumsum(cf, dims = 1)

    # Discounted NPV
    npv = sum(cf, dims=1)

    return Costs(capex, opex, salvage, cf, cumulative_npv, npv)
end

# Baseline cost
compute_baseline_cost(des::DistributedEnergySystem) = dropdims(sum(max.(0,  (isa(des.ld_E, Load) ? des.ld_E.power : 0. ) .+ (isa(des.ld_H, Load) ? des.ld_H.power ./ des.heater.η_E_H : 0.)) .* des.grid.cost_in .* des.parameters.Δh, dims=1), dims=1)
compute_baseline_cost(s::Int64, des::DistributedEnergySystem) = dropdims(sum(max.(0, (isa(des.ld_E, Load) ? des.ld_E.power[:,:,s] : 0. ) .+ (isa(des.ld_H, Load) ? des.ld_H.power[:,:,s] ./ des.heater.η_E_H : 0.)) .* des.grid.cost_in[:,:,s] .* des.parameters.Δh, dims=1), dims=1)
# Grid cost
compute_grid_cost(des::DistributedEnergySystem) = dropdims(sum((max.(0., des.grid.power_E) .* des.grid.cost_in .- min.(0., des.grid.power_E) .* des.grid.cost_out) .* des.parameters.Δh, dims=1), dims=1)
compute_grid_cost(s::Int64, des::DistributedEnergySystem) = dropdims(sum((max.(0., des.grid.power_E[:,:,s]) .* des.grid.cost_in[:,:,s] .- min.(0., des.grid.power_E[:,:,s]) .* des.grid.cost_out[:,:,s]) .* des.parameters.Δh, dims=1), dims=1)
compute_grid_cost(y::Int64, s::Int64, des::DistributedEnergySystem) = sum((max.(0., des.grid.power_E[:,y,s]) .* des.grid.cost_in[:,y,s] .- min.(0., des.grid.power_E[:,y,s]) .* des.grid.cost_out[:,y,s]) .* des.parameters.Δh)

# CAPEX
function compute_capex(des::DistributedEnergySystem, designer::AbstractDesigner)
    return designer.u.pv .* (isa(des.pv, Source) ? des.pv.cost : 0.) .+
            designer.u.liion .* (isa(des.liion, Liion) ? des.liion.cost : 0.) .+
            designer.u.tes .* (isa(des.tes, ThermalSto) ? des.tes.cost : 0.) .+
            designer.u.h2tank .* (isa(des.h2tank, H2Tank) ? des.h2tank.cost : 0.) .+
            designer.u.elyz .* (isa(des.elyz, Electrolyzer) ? des.elyz.cost : 0.) .+
            designer.u.fc .* (isa(des.fc, FuelCell) ? des.fc.cost : 0.)
end
# CAPEX for a given scenario s
function compute_capex(s::Int64, des::DistributedEnergySystem, designer::AbstractDesigner)
    return designer.u.pv[:,s] .* (isa(des.pv, Source) ? des.pv.cost[:,s] : 0.) .+
            designer.u.liion[:,s] .* (isa(des.liion, Liion) ? des.liion.cost[:,s] : 0.) .+
            designer.u.tes[:,s] .* (isa(des.tes, ThermalSto) ? des.tes.cost[:,s] : 0.) .+
            designer.u.h2tank[:,s] .* (isa(des.h2tank, H2Tank) ? des.h2tank.cost[:,s] : 0.) .+
            designer.u.elyz[:,s] .* (isa(des.elyz, Electrolyzer) ? des.elyz.cost[:,s] : 0.) .+
            designer.u.fc[:,s] .* (isa(des.fc, FuelCell) ? des.fc.cost[:,s] : 0.)
end
# Annualised CAPEX
function compute_annualised_capex(y::Int64, s::Int64, des::DistributedEnergySystem, designer::AbstractDesigner)
    # Preallocation
    capex = 0.

    if isa(des.pv, Source)
        Γ_pv = (des.parameters.τ * (des.parameters.τ + 1.) ^ des.pv.lifetime) / ((des.parameters.τ + 1.) ^ des.pv.lifetime - 1.)
        capex += Γ_pv * designer.u.pv[y,s] * des.pv.cost[y,s]
    end

    if isa(des.liion, Liion)
        Γ_liion = (des.parameters.τ * (des.parameters.τ + 1.) ^ des.liion.lifetime) / ((des.parameters.τ + 1.) ^ des.liion.lifetime - 1.)
        capex += Γ_liion * designer.u.liion[y,s] * des.liion.cost[y,s]
    end

    if isa(des.tes, ThermalSto)
        Γ_tes = (des.parameters.τ * (des.parameters.τ + 1.) ^ des.tes.lifetime) / ((des.parameters.τ + 1.) ^ des.tes.lifetime - 1.)
        capex += Γ_tes * designer.u.tes[y,s] * des.tes.cost[y,s]
    end

    if isa(des.h2tank, H2Tank)
        Γ_h2tank = (des.parameters.τ * (des.parameters.τ + 1.) ^ des.h2tank.lifetime) / ((des.parameters.τ + 1.) ^ des.h2tank.lifetime - 1.)
        capex += Γ_h2tank * designer.u.h2tank[y,s] * des.h2tank.cost[y,s]
    end

    if isa(des.elyz, Electrolyzer)
        Γ_elyz = (des.parameters.τ * (des.parameters.τ + 1.) ^ des.elyz.lifetime) / ((des.parameters.τ + 1.) ^ des.elyz.lifetime - 1.)
        capex += Γ_elyz * designer.u.elyz[y,s] * des.elyz.cost[y,s]
    end

    if isa(des.fc, FuelCell)
        Γ_fc = (des.parameters.τ * (des.parameters.τ + 1.) ^ des.fc.lifetime) / ((des.parameters.τ + 1.) ^ des.fc.lifetime - 1.)
        capex += Γ_fc * designer.u.fc[y,s] * des.fc.cost[y,s]
    end

    return capex
end

# Salvage value
function compute_salvage(des::DistributedEnergySystem)
    # Linear depreciation of components
    nh, ny = des.parameters.nh, des.parameters.ny
    salvage = zeros(des.parameters.ny, des.parameters.ns)

    salvage[ny, :] .= (isa(des.pv, Source) ? (des.pv.lifetime .- ny) ./ des.pv.lifetime .* des.pv.cost[ny, :] : 0.) .+
                       (isa(des.liion, Liion) ? des.liion.soh[nh, ny, :] .* des.liion.Erated[ny,:] .* des.liion.cost[ny, :] : 0.) .+
                       (isa(des.tes, ThermalSto) ? (des.tes.lifetime .- ny) ./ des.tes.lifetime .* des.tes.cost[ny, :] : 0.) .+
                       (isa(des.h2tank, H2Tank) ? (des.h2tank.lifetime .- ny) ./ des.h2tank.lifetime .* des.h2tank.cost[ny, :] : 0.) .+
                       (isa(des.elyz, Electrolyzer) ? des.elyz.soh[nh, ny, :] .* des.elyz.powerMax[ny,:] .* des.elyz.cost[ny, :] : 0.) .+
                       (isa(des.fc, FuelCell) ? des.fc.soh[nh, ny, :] .* des.fc.powerMax[ny,:] .* des.fc.cost[ny, :] : 0.)
    return salvage

end
# Salvage value for a given scenario s
function compute_salvage(s::Int64, des::DistributedEnergySystem)
    # Linear depreciation of components
    nh, ny = des.parameters.nh, des.parameters.ny
    salvage = zeros(des.parameters.ny)

    salvage[ny] = (isa(des.pv, Source) ? (des.pv.lifetime .- ny) ./ des.pv.lifetime .* des.pv.cost[ny, s] : 0.) +
                       (isa(des.liion, Liion) ? des.liion.soh[nh, ny, s] .* des.liion.Erated[ny,s] .* des.liion.cost[ny, s] : 0.) +
                       (isa(des.tes, ThermalSto) ? (des.tes.lifetime .- ny) ./ des.tes.lifetime .* des.tes.cost[ny, s] : 0.) +
                       (isa(des.h2tank, H2Tank) ? (des.h2tank.lifetime .- ny) ./ des.h2tank.lifetime .* des.h2tank.cost[ny, s] : 0.) +
                       (isa(des.elyz, Electrolyzer) ? des.elyz.soh[nh, ny, s] .* des.elyz.powerMax[ny,s] .* des.elyz.cost[ny, s] : 0.) +
                       (isa(des.fc, FuelCell) ? des.fc.soh[nh, ny, s] .* des.fc.powerMax[ny,s] .* des.fc.cost[ny, s] : 0.)
    return salvage
end

# Share of renewables
function compute_share(des::DistributedEnergySystem)
    # Total demand
    sum_ld_tot = dropdims(sum((isa(des.ld_E, Load) ? des.ld_E.power : 0. ) .+ (isa(des.ld_H, Load) ? des.ld_H.power ./ des.heater.η_E_H : 0.), dims=1), dims=1)

    # Share of renew.
    τ_share = 1. .- dropdims(sum(max.(0., des.grid.power_E), dims=1), dims=1) ./ sum_ld_tot

    return τ_share
end
# Share of renewables for a given scenario s
function compute_share(s::Int64, des::DistributedEnergySystem)
    # Total demand
    sum_ld_tot = dropdims(sum((isa(des.ld_E, Load) ? des.ld_E.power[:,:,s] : 0. ) .+ (isa(des.ld_H, Load) ? des.ld_H.power[:,:,s] ./ des.heater.η_E_H : 0.), dims=1), dims=1)

    # Share of renew.
    τ_share = 1. .- dropdims(sum(max.(0., des.grid.power_E[:,:,s]), dims=1), dims=1) ./ sum_ld_tot

    return τ_share
end
# Share of renewables for a given year y of a givn scenario s
function compute_share(y::Int64, s::Int64, des::DistributedEnergySystem)
    # Total demand
    sum_ld_tot = sum((isa(des.ld_E, Load) ? sum(des.ld_E.power[:,y,s]) : 0.) .+ (isa(des.ld_H, Load) ? sum(des.ld_H.power[:,y,s]) ./ des.heater.η_E_H : 0.))

    # Share of renew.
    τ_share = 1. .- sum(max.(0., des.grid.power_E[:,y,s])) ./ sum_ld_tot

    return τ_share
end

mutable struct LPSP
    lpsp_E::Union{Nothing, Float64, Array{Float64,1}, Array{Float64,2}}
    lpsp_H::Union{Nothing, Float64, Array{Float64,1}, Array{Float64,2}}
end

# LPSP
function LPSP(des::DistributedEnergySystem)

    # Elec.
    isa(des.ld_E, Load) ? ld_E = des.ld_E.power : ld_E = 0.
    isa(des.pv, Source) ? pv = des.pv.power_E : pv = 0.
    isa(des.heater, Heater) ? heater = des.heater.power_E : heater = 0.
    isa(des.liion, Liion) ? liion = des.liion.power_E : liion = 0.
    isa(des.elyz, Electrolyzer) ? elyz = des.elyz.power_E : elyz = 0.
    isa(des.fc, FuelCell) ? fc = des.fc.power_E : fc = 0.
    isa(des.grid, Grid) ? grid = des.grid.power_E : grid = 0.

    isa(des.ld_E, Load) ? lpsp_E = dropdims(sum(max.(0., ld_E .- pv .- liion .- elyz .- fc .- heater .- grid), dims=1) ./ sum(ld_E, dims=1), dims=1) : lpsp_E = nothing

    # Heat
    isa(des.ld_H, Load) ? ld_H = des.ld_H.power : ld_H = 0.
    isa(des.heater, Heater) ? heater = des.heater.power_H : heater = 0.
    isa(des.tes, ThermalSto) ? tes = des.tes.power_H : tes = 0.
    isa(des.elyz, Electrolyzer) ? elyz = des.elyz.power_H : elyz = 0.
    isa(des.fc, FuelCell) ? fc = des.fc.power_H : fc = 0.

    isa(des.ld_H, Load) ? lpsp_H = dropdims(sum(max.(0., ld_H .- heater .- fc .- elyz .- tes), dims=1) ./ sum(ld_H, dims=1), dims=1) : lpsp_H = nothing

    return LPSP(lpsp_E, lpsp_H)
end
# LPSP for a given scenario s
function LPSP(s::Int64, des::DistributedEnergySystem)

    # Elec.
    isa(des.ld_E, Load) ? ld_E = des.ld_E.power[:,:,s] : ld_E = 0.
    isa(des.pv, Source) ? pv = des.pv.power_E[:,:,s] : pv = 0.
    isa(des.heater, Heater) ? heater = des.heater.power_E[:,:,s] : heater = 0.
    isa(des.liion, Liion) ? liion = des.liion.power_E[:,:,s] : liion = 0.
    isa(des.elyz, Electrolyzer) ? elyz = des.elyz.power_E[:,:,s] : elyz = 0.
    isa(des.fc, FuelCell) ? fc = des.fc.power_E[:,:,s] : fc = 0.
    isa(des.grid, Grid) ? grid = des.grid.power_E[:,:,s] : grid = 0.

    isa(des.ld_E, Load) ? lpsp_E = dropdims(sum(max.(0., ld_E .- pv .- liion .- elyz .- fc .- heater .- grid), dims=1) ./ sum(ld_E, dims=1), dims=1) : lpsp_E = nothing

    # Heat
    isa(des.ld_H, Load) ? ld_H = des.ld_H.power[:,:,s] : ld_H = 0.
    isa(des.heater, Heater) ? heater = des.heater.power_H[:,:,s] : heater = 0.
    isa(des.tes, ThermalSto) ? tes = des.tes.power_H[:,:,s] : tes = 0.
    isa(des.elyz, Electrolyzer) ? elyz = des.elyz.power_H[:,:,s] : elyz = 0.
    isa(des.fc, FuelCell) ? fc = des.fc.power_H[:,:,s] : fc = 0.

    isa(des.ld_H, Load) ? lpsp_H = dropdims(sum(max.(0., ld_H .- heater .- fc .- elyz .- tes), dims=1) ./ sum(ld_H, dims=1), dims=1) : lpsp_H = nothing

    return LPSP(lpsp_E, lpsp_H)
end
# LPSP for a given scenario s and year y
function LPSP(y::Int64, s::Int64, des::DistributedEnergySystem)

    # Elec.
    isa(des.ld_E, Load) ? ld_E = des.ld_E.power[:,y,s] : ld_E = 0.
    isa(des.pv, Source) ? pv = des.pv.power_E[:,y,s] : pv = 0.
    isa(des.heater, Heater) ? heater = des.heater.power_E[:,y,s] : heater = 0.
    isa(des.liion, Liion) ? liion = des.liion.power_E[:,y,s] : liion = 0.
    isa(des.elyz, Electrolyzer) ? elyz = des.elyz.power_E[:,y,s] : elyz = 0.
    isa(des.fc, FuelCell) ? fc = des.fc.power_E[:,y,s] : fc = 0.
    isa(des.grid, Grid) ? grid = des.grid.power_E[:,y,s] : grid = 0.

    isa(des.ld_E, Load) ? lpsp_E = sum(max.(0., ld_E .- pv .- liion .- elyz .- fc .- heater .- grid)) ./ sum(ld_E) : lpsp_E = nothing

    # Heat
    isa(des.ld_H, Load) ? ld_H = des.ld_H.power[:,y,s] : ld_H = 0.
    isa(des.heater, Heater) ? heater = des.heater.power_H[:,y,s] : heater = 0.
    isa(des.tes, ThermalSto) ? tes = des.tes.power_H[:,y,s] : tes = 0.
    isa(des.elyz, Electrolyzer) ? elyz = des.elyz.power_H[:,y,s] : elyz = 0.
    isa(des.fc, FuelCell) ? fc = des.fc.power_H[:,y,s] : fc = 0.

    isa(des.ld_H, Load) ? lpsp_H = sum(max.(0., ld_H .- heater .- fc .- elyz .- tes)) ./ sum(ld_H) : lpsp_H = nothing

    return LPSP(lpsp_E, lpsp_H)
end

mutable struct Metrics
    costs::Costs
    τ_share::Union{Array{Float64,1}, Array{Float64,2}}
    lpsp::LPSP
end

# Compute indicators
function Metrics(des::DistributedEnergySystem, designer::AbstractDesigner)
    # Econmics
    costs = Costs(des, designer)

    # Share of renewables
    τ_share = compute_share(des)

    # LPSP
    lpsp = LPSP(des)

    return Metrics(costs, τ_share, lpsp)
end
# Compute indicators for a given scenario s
function Metrics(s::Int64, des::DistributedEnergySystem, designer::AbstractDesigner)
    # Econmics
    costs = Costs(s, des, designer)

    # Share of renewables
    τ_share = compute_share(s, des)

    # LPSP
    lpsp = LPSP(s, des)

    return Metrics(costs, τ_share, lpsp)
end
