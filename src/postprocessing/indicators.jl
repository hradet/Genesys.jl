#=
    This file includes all the funtions needed to compute the techno-economic
    indicators
 =#

### Economics
 mutable struct Costs
      capex::Array{Float64,2}
      opex::Array{Float64,2}
      cf::Array{Float64,2}
      cumulative_npv::Array{Float64,2}
      npv::Array{Float64,1}
 end

function compute_economics(des::DistributedEnergySystem, designer::AbstractDesigner)

     # Discount factor
     γ = 1. ./ (1. + des.parameters.τ) .^ range(1, length = des.parameters.ny, step = des.parameters.Δy)

    # Discounted capex
    capex = γ .* compute_capex(des, designer)

    # Discounted opex
    opex = γ .* compute_opex(des)

    # Discounted cash flow
    cf = - capex .+ opex

    # Discounted NPV each year
    cumulative_npv = cumsum(cf, dims = 1)

    # Discounted NPV
    npv = dropdims(sum(cf, dims=1), dims=1)

    return Costs(capex, opex, cf, cumulative_npv, npv)
end
function compute_opex(des::DistributedEnergySystem)
    # Total load
    ld_tot = zeros(des.parameters.nh, des.parameters.ny, des.parameters.ns)
    isa(des.ld_E, Load) ? ld_tot .+= des.ld_E.power : nothing
    isa(des.ld_H, Load) ? ld_tot .+= des.ld_H.power ./ des.heater.η_E_H : nothing

    # Reference case when all the electricity is purchased from the grid
    ref = max.(0, ld_tot)
    saving = ((ref .- max.(0., des.grid.power_E)) .* des.grid.C_grid_in - min.(0., des.grid.power_E) .* des.grid.C_grid_out) .* des.parameters.Δh

    # opex
    opex = dropdims(sum(saving, dims=1), dims=1)

    return opex
end
function compute_capex(des::DistributedEnergySystem, designer::AbstractDesigner)
    # Preallocation
    capex = zeros(des.parameters.ny, des.parameters.ns)

    isa(des.pv, Source) ? capex .+= designer.u.pv .* des.pv.C_pv : nothing
    isa(des.liion, Liion) ? capex .+= designer.u.liion .* des.liion.C_liion : nothing
    isa(des.tes, ThermalSto) ? capex .+= designer.u.tes .* des.tes.C_tes : nothing
    isa(des.h2tank, H2Tank) ? capex .+= designer.u.h2tank .* des.h2tank.C_tank : nothing
    isa(des.elyz, Electrolyzer) ? capex .+= designer.u.elyz .* des.elyz.C_elyz : nothing
    isa(des.fc, FuelCell) ? capex .+= designer.u.fc .* des.fc.C_fc : nothing

    return capex
end

### Techno ###
mutable struct Indicators
     τ_share::Array{Float64,2}
     τ_autoconso
end

function compute_tech_indicators(des::DistributedEnergySystem)
    # Self-sufficiency
    τ_share = dropdims(1. .- sum(max.(0., des.grid.power_E), dims=1) ./ sum(des.ld_E.power, dims=1), dims=1)
    # Self-consumption
    # TODO

    return Indicators(τ_share, nothing)
end
