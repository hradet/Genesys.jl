#=
    This file includes all the funtions needed to compute the economic
    indicators
 =#
 
 # Cost of operation
function compute_opex(ld::Load, grid::Grid, Δh)

    # Reference case when all the electricity is purchased from the grid
    ref = max.(0,ld.power_E)
    saving = (ref .- grid.power_E) .* Δh .* grid.C_grid_in

    # opex
    opex = dropdims(sum(saving, dims=1), dims=1)

    return opex
end
# Cost of investment
function compute_capex(pv::Source, liion::Liion, designer::AbstractDesigner)

    #capex
    capex = designer.u.u_pv .* pv.C_pv .+ designer.u.u_liion .* liion.C_liion

    return capex
end
# Total cost
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
