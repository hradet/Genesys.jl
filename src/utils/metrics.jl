#=
    This file includes all the funtions needed to compute the techno-economic
    indicators
 =#


 mutable struct Costs{T <: Array{Float64}}
      capex::T
      opex::T
      salvage::T
      cf::T
      cumulative_npv::T
      npv::T
 end

 # Compute costs
Costs(des::DistributedEnergySystem, designer::AbstractDesigner) = Costs(1:des.parameters.ns, des, designer)
# Compute costs for a given scenario s
function Costs(s::Union{Int64, UnitRange{Int64}}, des::DistributedEnergySystem, designer::AbstractDesigner)

     # Discount factor
     γ = repeat(1. ./ (1. + des.parameters.τ) .^ range(0, length = des.parameters.ny, step = des.parameters.Δy), 1, length(s))

    # Discounted capex
    capexx = γ .* capex(s, des, designer)

    # Discounted opex
    opex = γ .* (baseline_cost(s, des) .- grid_cost(s, des))

    # Discounted salvage value
    salvage = γ .* salvage_value(s, des)

    # Discounted cash flow
    cf = - capexx .+ opex .+ salvage

    # Discounted NPV each year
    cumulative_npv = cumsum(cf, dims = 1)

    # Discounted NPV
    npv = sum(cf, dims=1)

    return Costs(capexx, opex, salvage, cf, cumulative_npv, npv)
end

# Baseline cost
baseline_cost(des::DistributedEnergySystem) = baseline_cost(1:des.parameters.ns, des)
baseline_cost(s::Union{Int64, UnitRange{Int64}}, des::DistributedEnergySystem) = sum(max.(0, (isa(des.ld_E, Load) ? des.ld_E.power[:,:,s] : 0. ) .+ (isa(des.ld_H, Load) ? des.ld_H.power[:,:,s] ./ des.heater.η_E_H : 0.)) .* des.grid.cost_in[:,:,s] .* des.parameters.Δh, dims=1)[1,:,:]
# Grid cost
grid_cost(des::DistributedEnergySystem) = grid_cost(1:des.parameters.ns, des)
grid_cost(s::Union{Int64, UnitRange{Int64}}, des::DistributedEnergySystem) = grid_cost(1:des.parameters.ny, s, des)
grid_cost(y::Union{Int64, UnitRange{Int64}}, s::Union{Int64, UnitRange{Int64}}, des::DistributedEnergySystem) = sum((max.(0., des.grid.power_E[:,y,s]) .* des.grid.cost_in[:,y,s] .- min.(0., des.grid.power_E[:,y,s]) .* des.grid.cost_out[:,y,s]) .* des.parameters.Δh, dims=1)[1,:,:]

# CAPEX
capex(des::DistributedEnergySystem, designer::AbstractDesigner) = capex(1:des.parameters.ns, des, designer)
# CAPEX for a given scenario s
function capex(s::Union{Int64, UnitRange{Int64}}, des::DistributedEnergySystem, designer::AbstractDesigner)
    return designer.u.pv[:,s] .* (isa(des.pv, Source) ? des.pv.cost[:,s] : 0.) .+
            designer.u.liion[:,s] .* (isa(des.liion, Liion) ? des.liion.cost[:,s] : 0.) .+
            designer.u.tes[:,s] .* (isa(des.tes, ThermalSto) ? des.tes.cost[:,s] : 0.) .+
            designer.u.h2tank[:,s] .* (isa(des.h2tank, H2Tank) ? des.h2tank.cost[:,s] : 0.) .+
            designer.u.elyz[:,s] .* (isa(des.elyz, Electrolyzer) ? des.elyz.cost[:,s] : 0.) .+
            designer.u.fc[:,s] .* (isa(des.fc, FuelCell) ? des.fc.cost[:,s] : 0.)
end
# Annualised CAPEX
function annualised_capex(y::Union{Int64, UnitRange{Int64}}, s::Union{Int64, UnitRange{Int64}}, des::DistributedEnergySystem, designer::AbstractDesigner)
    # Preallocation
    capex = 0.

    if isa(des.pv, Source)
        Γ_pv = (des.parameters.τ * (des.parameters.τ + 1.) ^ des.pv.lifetime) / ((des.parameters.τ + 1.) ^ des.pv.lifetime - 1.)
        capex = capex .+ Γ_pv .* designer.u.pv[y,s] .* des.pv.cost[y,s]
    end

    if isa(des.liion, Liion)
        Γ_liion = (des.parameters.τ * (des.parameters.τ + 1.) ^ des.liion.lifetime) / ((des.parameters.τ + 1.) ^ des.liion.lifetime - 1.)
        capex = capex .+ Γ_liion .* designer.u.liion[y,s] .* des.liion.cost[y,s]
    end

    if isa(des.tes, ThermalSto)
        Γ_tes = (des.parameters.τ * (des.parameters.τ + 1.) ^ des.tes.lifetime) / ((des.parameters.τ + 1.) ^ des.tes.lifetime - 1.)
        capex = capex .+ Γ_tes .* designer.u.tes[y,s] .* des.tes.cost[y,s]
    end

    if isa(des.h2tank, H2Tank)
        Γ_h2tank = (des.parameters.τ * (des.parameters.τ + 1.) ^ des.h2tank.lifetime) / ((des.parameters.τ + 1.) ^ des.h2tank.lifetime - 1.)
        capex = capex .+ Γ_h2tank .* designer.u.h2tank[y,s] .* des.h2tank.cost[y,s]
    end

    if isa(des.elyz, Electrolyzer)
        Γ_elyz = (des.parameters.τ * (des.parameters.τ + 1.) ^ des.elyz.lifetime) / ((des.parameters.τ + 1.) ^ des.elyz.lifetime - 1.)
        capex = capex .+ Γ_elyz .* designer.u.elyz[y,s] .* des.elyz.cost[y,s]
    end

    if isa(des.fc, FuelCell)
        Γ_fc = (des.parameters.τ * (des.parameters.τ + 1.) ^ des.fc.lifetime) / ((des.parameters.τ + 1.) ^ des.fc.lifetime - 1.)
        capex = capex .+ Γ_fc .* designer.u.fc[y,s] .* des.fc.cost[y,s]
    end

    return capex
end

# Salvage value
salvage_value(des::DistributedEnergySystem) = salvage_value(1:des.parameters.ns, des)
# Salvage value for a given scenario s
function salvage_value(s::Union{Int64, UnitRange{Int64}}, des::DistributedEnergySystem)
    # Linear depreciation of components
    nh, ny = des.parameters.nh, des.parameters.ny
    salvage = zeros(des.parameters.ny, length(s))

    salvage[ny,:] .= (isa(des.pv, Source) ? (des.pv.lifetime .- ny) ./ des.pv.lifetime .* des.pv.cost[ny, s] : 0.) .+
                       (isa(des.liion, Liion) ? des.liion.soh[nh, ny, s] .* des.liion.Erated[ny,s] .* des.liion.cost[ny, s] : 0.) .+
                       (isa(des.tes, ThermalSto) ? (des.tes.lifetime .- ny) ./ des.tes.lifetime .* des.tes.cost[ny, s] : 0.) .+
                       (isa(des.h2tank, H2Tank) ? (des.h2tank.lifetime .- ny) ./ des.h2tank.lifetime .* des.h2tank.cost[ny, s] : 0.) .+
                       (isa(des.elyz, Electrolyzer) ? des.elyz.soh[nh, ny, s] .* des.elyz.powerMax[ny,s] .* des.elyz.cost[ny, s] : 0.) .+
                       (isa(des.fc, FuelCell) ? des.fc.soh[nh, ny, s] .* des.fc.powerMax[ny,s] .* des.fc.cost[ny, s] : 0.)
    return salvage
end

# Share of renewables
renewable_share(des::DistributedEnergySystem) = renewable_share(1:des.parameters.ns, des)
# Share of renewables for a given scenario s
renewable_share(s::Union{Int64, UnitRange{Int64}}, des::DistributedEnergySystem) = renewable_share(1:des.parameters.ny, s, des)
# Share of renewables for a given year y of a givn scenario s
function renewable_share(y::Union{Int64, UnitRange{Int64}}, s::Union{Int64, UnitRange{Int64}}, des::DistributedEnergySystem)
    return (1. .- sum(max.(0., des.grid.power_E[:,y,s]), dims = 1) ./ ((isa(des.ld_E, Load) ? sum(des.ld_E.power[:,y,s], dims = 1) : 0.) .+ (isa(des.ld_H, Load) ? sum(des.ld_H.power[:,y,s], dims = 1) ./ des.heater.η_E_H : 0.)))[1,:,:]
end

# LPSP
mutable struct LPSP{T}
    lpsp_E::Union{Nothing, T}
    lpsp_H::Union{Nothing, T}
end

LPSP(des::DistributedEnergySystem) = LPSP(1:des.parameters.ns, des)
# LPSP for a given scenario s
LPSP(s::Union{Int64, UnitRange{Int64}}, des::DistributedEnergySystem) = LPSP(1:des.parameters.ny, s, des)
# LPSP for a given scenario s and year y
function LPSP(y::Union{Int64, UnitRange{Int64}}, s::Union{Int64, UnitRange{Int64}}, des::DistributedEnergySystem)

    # Elec.
    isa(des.ld_E, Load) ? ld_E = des.ld_E.power[:,y,s] : ld_E = 0.
    isa(des.pv, Source) ? pv = des.pv.power_E[:,y,s] : pv = 0.
    isa(des.heater, Heater) ? heater = des.heater.power_E[:,y,s] : heater = 0.
    isa(des.liion, Liion) ? liion = des.liion.power_E[:,y,s] : liion = 0.
    isa(des.elyz, Electrolyzer) ? elyz = des.elyz.power_E[:,y,s] : elyz = 0.
    isa(des.fc, FuelCell) ? fc = des.fc.power_E[:,y,s] : fc = 0.
    isa(des.grid, Grid) ? grid = des.grid.power_E[:,y,s] : grid = 0.

    isa(des.ld_E, Load) ? lpsp_E = sum(max.(0., ld_E .- pv .- liion .- elyz .- fc .- heater .- grid), dims=1)[1,:,:] ./ sum(ld_E, dims=1)[1,:,:] : lpsp_E = nothing

    # Heat
    isa(des.ld_H, Load) ? ld_H = des.ld_H.power[:,y,s] : ld_H = 0.
    isa(des.heater, Heater) ? heater = des.heater.power_H[:,y,s] : heater = 0.
    isa(des.tes, ThermalSto) ? tes = des.tes.power_H[:,y,s] : tes = 0.
    isa(des.elyz, Electrolyzer) ? elyz = des.elyz.power_H[:,y,s] : elyz = 0.
    isa(des.fc, FuelCell) ? fc = des.fc.power_H[:,y,s] : fc = 0.

    isa(des.ld_H, Load) ? lpsp_H = sum(max.(0., ld_H .- heater .- fc .- elyz .- tes), dims=1)[1,:,:] ./ sum(ld_H, dims=1)[1,:,:] : lpsp_H = nothing

    return LPSP(lpsp_E, lpsp_H)
end

mutable struct Metrics{T}
    costs::Costs{T}
    renewable_share::T
    lpsp::LPSP{T}
end

# Compute indicators
Metrics(des::DistributedEnergySystem, designer::AbstractDesigner) = Metrics(1:des.parameters.ns, des, designer)
# Compute indicators for a given scenario s
function Metrics(s::Union{Int64, UnitRange{Int64}}, des::DistributedEnergySystem, designer::AbstractDesigner)
    # Econmics
    costs = Costs(s, des, designer)

    # Share of renewables
    share = renewable_share(s, des)

    # LPSP
    lpsp = LPSP(s, des)

    return Metrics(costs, share, lpsp)
end
