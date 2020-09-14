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

function compute_economics(des::DES)

     # Discount factor
     γ = 1. ./ (1. + des.parameters.τ) .^ (1:des.parameters.Y)

    # Investment cost
    capex = compute_capex(des)

    # Operation cost
    opex = compute_opex(des)

    # Cash flow
    cf = γ .* (- capex .+ opex)

    # Discounted NPV each year
    cumulative_npv = cumsum(cf, dims = 1)

    # Discounted NPV
    npv = dropdims(sum(cf, dims=1), dims=1)

    return Costs(capex, opex, cf, cumulative_npv, npv)
end
function compute_opex(des::DES)
    # Parameters
    nh = size(des.ld_E.power, 1)
    ny = size(des.ld_E.power, 2)
    ns = size(des.ld_E.power, 3)

    isa(des.ld_H, Load) ? ld_H_to_E = des.ld_H.power ./ des.heater.η_E_H : ld_H_to_E = zeros(nh, ny, ns)

    # Reference case when all the electricity is purchased from the grid
    ref = max.(0, des.ld_E.power .+ ld_H_to_E)
    saving = ((ref .- max.(0., des.grid.power_E)) .* des.grid.C_grid_in - min.(0., des.grid.power_E) .* des.grid.C_grid_out) .* des.parameters.Δh

    # opex
    opex = dropdims(sum(saving, dims=1), dims=1)

    return opex
end
function compute_capex(des::DES)
    # Parameters
    ny = size(des.ld_E.power, 2)
    ns = size(des.ld_E.power, 3)
    # Preallocation
    capex = zeros(ny, ns)

    isa(des.pv, Source) ? capex .+= des.designer.u.pv .* des.pv.C_pv : nothing
    isa(des.liion, Liion) ? capex .+= des.designer.u.liion .* des.liion.C_liion : nothing
    isa(des.tes, ThermalSto) ? capex .+= des.designer.u.tes .* des.tes.C_tes : nothing
    isa(des.h2tank, H2Tank) ? capex .+= des.designer.u.h2tank .* des.h2tank.C_tank : nothing
    isa(des.elyz, Electrolyzer) ? capex .+= des.designer.u.elyz .* des.elyz.C_elyz : nothing
    isa(des.fc, FuelCell) ? capex .+= des.designer.u.fc .* des.fc.C_fc : nothing

    return capex
end

### Techno ###
mutable struct Indicators
     τ_self::Array{Float64,2}
     τ_autoconso
end

function compute_tech_indicators(des::DES)
    # Self-sufficiency
    τ_self = dropdims(1. .- sum(max.(0., des.grid.power_E), dims=1) ./ sum(des.ld_E.power, dims=1), dims=1)
    # Self-consumption
    # TODO

    return Indicators(τ_self, nothing)
end
