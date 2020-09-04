#=
    This file includes all the funtions needed to compute the economic
    indicators
 =#

 mutable struct Costs
      capex::Array{Float64,2}
      opex::Array{Float64,2}
      cf::Array{Float64,2}
      cumulative_npv::Array{Float64,2}
      npv::Array{Float64,1}
 end

# Simple
function compute_economics(ld::Load, pv::Source, liion::Liion, designer::AbstractDesigner,
     grid::Grid, parameters::NamedTuple)

    # Discount factor
    γ = 1. ./ (1. + parameters.τ) .^ (1:parameters.Y)

    # Investment cost
    capex = compute_capex(pv, liion, designer)

    # Operation cost
    opex = compute_opex(ld, grid, parameters.Δh)

    # Cash flow
    cf = γ .* (- capex .+ opex)

    # Discounted NPV each year
    cumulative_npv = cumsum(cf, dims = 1)

    # Discounted NPV
    npv = dropdims(sum(cf, dims=1), dims=1)

    return Costs(capex, opex, cf, cumulative_npv, npv)
end
function compute_opex(ld::Load, grid::Grid, Δh)

    # Reference case when all the electricity is purchased from the grid
    ref = max.(0,ld.power_E)
    saving = (ref .- grid.power_E) .* Δh .* grid.C_grid_in

    # opex
    opex = dropdims(sum(saving, dims=1), dims=1)

    return opex
end
function compute_capex(pv::Source, liion::Liion, designer::AbstractDesigner)

    #capex
    capex = designer.u.u_pv .* pv.C_pv .+ designer.u.u_liion .* liion.C_liion

    return capex
end

# Multi-energy
function compute_economics(ld::Load, pv::Source, liion::Liion, h2tank::H2Tank,
     elyz::Electrolyzer, fc::FuelCell, tes::ThermalSto, heater::Heater,
     designer::AbstractDesigner, grid::Grid, parameters::NamedTuple)

     # Discount factor
     γ = 1. ./ (1. + parameters.τ) .^ (1:parameters.Y)

    # Investment cost
    capex = compute_capex(pv, liion, h2tank, elyz, fc, tes, designer)

    # Operation cost
    opex = compute_opex(ld, heater, grid, parameters.Δh)

    # Cash flow
    cf = γ .* (- capex .+ opex)

    # Discounted NPV each year
    cumulative_npv = cumsum(cf, dims = 1)

    # Discounted NPV
    npv = dropdims(sum(cf, dims=1), dims=1)

    return Costs(capex, opex, cf, cumulative_npv, npv)
end

function compute_opex(ld::Load, heater::Heater, grid::Grid, Δh)

    # Reference case when all the electricity is purchased from the grid
    ref = max.(0,ld.power_E .+ ld.power_H ./ heater.η_E_H)
    saving = (ref .- grid.power_E) .* Δh .* grid.C_grid_in

    # opex
    opex = dropdims(sum(saving, dims=1), dims=1)

    return opex
end
function compute_capex(pv::Source, liion::Liion, h2tank::H2Tank, elyz::Electrolyzer,
     fc::FuelCell, tes::ThermalSto, designer::AbstractDesigner)

    #capex
    capex = designer.u.u_pv .* pv.C_pv .+ designer.u.u_liion .* liion.C_liion .+
    designer.u.u_tank .* h2tank.C_tank .+ designer.u.u_elyz .* elyz.C_elyz .+
    designer.u.u_fc .* fc.C_fc .+ designer.u.u_tes .* tes.C_tes

    return capex
end
