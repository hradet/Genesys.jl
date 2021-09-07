#=
    This file includes all the funtions needed to compute the techno-economic
    indicators
 =#


 mutable struct NPV{T <: Array{Float64}}
      capex::T
      opex::T
      salvage::T
      cf::T
      cumulative_npv::T
      total::T
 end

 # Compute costs
NPV(mg::Microgrid, designer::AbstractDesigner) = NPV(1:mg.parameters.ns, mg, designer)
# Compute costs for a given scenario s
function NPV(s::Union{Int64, UnitRange{Int64}}, mg::Microgrid, designer::AbstractDesigner)

    # Discount factor
    γ = repeat(1. ./ (1. + mg.parameters.τ) .^ range(0, length = mg.parameters.ny, step = mg.parameters.Δy), 1, length(s))

    # Discounted capex
    capexx = γ .* capex(s, mg, designer)

    # Discounted opex
    opex = γ .* dropdims(baseline_cost(s, mg) .-  grid_cost(s, mg), dims=1)

    # Discounted salvage value
    salvage = γ .*  salvage_value(s, mg)

    # Discounted cash flow
    cf = - capexx .+ opex .+ salvage

    # Discounted NPV each year
    cumulative = cumsum(cf, dims = 1)

    # Discounted NPV
    total = sum(cf, dims=1)

    return NPV(capexx, opex, salvage, cf, cumulative, total)
end

# Baseline cost
baseline_cost(mg::Microgrid) = baseline_cost(1:mg.parameters.ns, mg)
baseline_cost(s::Union{Int64, UnitRange{Int64}}, mg::Microgrid) = baseline_cost(1:mg.parameters.ny, s, mg)
function baseline_cost(y::Union{Int64, UnitRange{Int64}}, s::Union{Int64, UnitRange{Int64}}, mg::Microgrid)
    elec, heat = 0., 0.
    if !isempty(mg.demands)
        for a in mg.demands
            if a.carrier isa Electricity
                elec = elec .+ sum(a.carrier.out[:,y,s] .* mg.grids[1].cost_in[:,y,s] * mg.parameters.Δh, dims = 1)
            elseif a.carrier isa Heat
                heat = heat .+ sum(a.carrier.out[:,y,s] / mg.converters[1].η_E_H .* mg.grids[1].cost_in[:,y,s] * mg.parameters.Δh, dims = 1)
            end
        end
    else
        println("Baseline cost does not exist!")
    end
    return elec .+ heat
end
# Grid cost
grid_cost(mg::Microgrid) = grid_cost(1:mg.parameters.ns, mg)
grid_cost(s::Union{Int64, UnitRange{Int64}}, mg::Microgrid) = grid_cost(1:mg.parameters.ny, s, mg)
grid_cost(y::Union{Int64, UnitRange{Int64}}, s::Union{Int64, UnitRange{Int64}}, mg::Microgrid) = sum(mg.grids[1].carrier.in[:,y,s] .* mg.grids[1].cost_in[:,y,s] .- mg.grids[1].carrier.out[:,y,s] .* mg.grids[1].cost_out[:,y,s], dims = 1) * mg.parameters.Δh

# CAPEX
capex(mg::Microgrid, designer::AbstractDesigner) = capex(1:mg.parameters.ns, mg, designer)
# CAPEX for a given scenario s
function capex(s::Union{Int64, UnitRange{Int64}}, mg::Microgrid, designer::AbstractDesigner)
    capex = 0.
    # Generations
    for (k, a) in enumerate(mg.generations)
        capex = capex .+ designer.decisions.generations[k][:,s] .* a.cost[:,s]
    end
    # Storages
    for (k, a) in enumerate(mg.storages)
        capex = capex .+ designer.decisions.storages[k][:,s] .* a.cost[:,s]
    end
    # Converters
    for (k, a) in enumerate(mg.converters)
        capex = capex .+ designer.decisions.converters[k][:,s] .* a.cost[:,s]
    end
    return capex
end

# Salvage value
salvage_value(mg::Microgrid) = salvage_value(1:mg.parameters.ns, mg)
# Salvage value for a given scenario s
function salvage_value(s::Union{Int64, UnitRange{Int64}}, mg::Microgrid)
    # TODO salvage as a function of SoH
    # Linear depreciation of components
    nh, ny = mg.parameters.nh, mg.parameters.ny
    salvage = zeros(mg.parameters.ny, length(s))
    salvage[ny,:] .= 0.
    # Generations
    for a in mg.generations
        salvage[ny,:] = salvage[ny,:] .+ (a.lifetime .- ny) ./ a.lifetime .* a.cost[ny, s]
    end
    # Storages
    for a in mg.storages
        salvage[ny,:] = salvage[ny,:] .+ (a.lifetime .- ny) ./ a.lifetime .* a.cost[ny, s]
    end
    # Converters
    for a in mg.converters
        salvage[ny,:] = salvage[ny,:] .+ (a.lifetime .- ny) ./ a.lifetime .* a.cost[ny, s]
    end
    return salvage
end

mutable struct EAC{T <: Array{Float64}}
     capex::T
     opex::T
     total::T
end

function EAC(y::Int64, s::Union{Int64, UnitRange{Int64}}, mg::Microgrid, designer::AbstractDesigner)
    # Annualised capex
    capex = Genesys.annualised_capex(y:y, s, mg, designer)
    # opex
    opex = grid_cost(y+1, s, mg)

    return EAC(capex, opex, capex .+ opex)
end
# Annualised CAPEX
function annualised_capex(y::Union{Int64, UnitRange{Int64}}, s::Union{Int64, UnitRange{Int64}}, mg::Microgrid, designer::AbstractDesigner)
    # Preallocation
    capex = 0.
    # Generations
    for (k, a) in enumerate(mg.generations)
        Γ = (mg.parameters.τ * (mg.parameters.τ + 1.) ^ a.lifetime) / ((mg.parameters.τ + 1.) ^ a.lifetime - 1.)
        capex = capex .+ Γ .* designer.decisions.generations[k][y,s] .* a.cost[y,s]
    end
    # Storages
    for (k, a) in enumerate(mg.storages)
        Γ = (mg.parameters.τ * (mg.parameters.τ + 1.) ^ a.lifetime) / ((mg.parameters.τ + 1.) ^ a.lifetime - 1.)
        capex = capex .+ Γ .* designer.decisions.storages[k][y,s] .* a.cost[y,s]
    end
    # Converters
    for (k, a) in enumerate(mg.converters)
        Γ = (mg.parameters.τ * (mg.parameters.τ + 1.) ^ a.lifetime) / ((mg.parameters.τ + 1.) ^ a.lifetime - 1.)
        capex = capex .+ Γ .* designer.decisions.converters[k][y,s] .* a.cost[y,s]
    end
    return capex
end

# Share of renewables
renewable_share(mg::Microgrid) = renewable_share(1:mg.parameters.ns, mg)
# Share of renewables for a given scenario s
renewable_share(s::Union{Int64, UnitRange{Int64}}, mg::Microgrid) = renewable_share(1:mg.parameters.ny, s, mg)
# Share of renewables for a given year y of a givn scenario s
function renewable_share(y::Union{Int64, UnitRange{Int64}}, s::Union{Int64, UnitRange{Int64}}, mg::Microgrid)
    # TODO to be changed if there is not grid...
    elec, heat = 0., 0.
    grid = sum(mg.grids[1].carrier.in[:,y,s], dims = 1)
    if !isempty(mg.demands)
        for a in mg.demands
            if a.carrier isa Electricity
                elec = elec .+ sum(a.carrier.out[:,y,s], dims = 1)
            elseif a.carrier isa Heat
                heat = heat .+ sum(a.carrier.out[:,y,s] / mg.converters[1].η_E_H , dims = 1)
            end
        end
    end
    return dropdims(1. .- grid ./ (elec .+ heat), dims=1)
end

# LPSP
mutable struct LPSP{T}
    elec::Union{Nothing, T}
    heat::Union{Nothing, T}
    hydrogen::Union{Nothing, T}
end

LPSP(mg::Microgrid) = LPSP(1:mg.parameters.ns, mg)
# LPSP for a given scenario s
LPSP(s::Union{Int64, UnitRange{Int64}}, mg::Microgrid) = LPSP(1:mg.parameters.ny, s, mg)
# LPSP for a given scenario s and year y
function LPSP(y::Union{Int64, UnitRange{Int64}}, s::Union{Int64, UnitRange{Int64}}, mg::Microgrid)
    # Initialization
    elec, heat, hydrogen = zeros(length(y), length(s)), zeros(length(y), length(s)), zeros(length(y), length(s))
    # Computation
    for a in mg.demands, ss in s, yy in y
        if a.carrier isa Electricity
            elec[yy,ss] = sum(max(0., power_balance(hh, yy, ss, mg, typeof(Electricity()))) for hh in 1:mg.parameters.nh) / sum(a.carrier.out[hh, yy, ss] for hh in 1:mg.parameters.nh)
            for aa in mg.grids
                if aa.carrier isa Electricity
                    elec[yy,ss] = elec[yy,ss] - sum(aa.carrier.in[hh, yy, ss] for hh in 1:mg.parameters.nh) / sum(a.carrier.out[hh, yy, ss] for hh in 1:mg.parameters.nh)
                end
            end
        elseif a.carrier isa Heat
            heat[yy,ss] = sum(max(0., power_balance(hh, yy, ss, mg, typeof(Heat()))) for hh in 1:mg.parameters.nh) / sum(a.carrier.out[hh, yy, ss] for hh in 1:mg.parameters.nh)
            for aa in mg.grids
                if aa.carrier isa Heat
                    heat[yy,ss] = heat[yy,ss] - sum(aa.carrier.in[hh, yy, ss] for hh in 1:mg.parameters.nh) / sum(a.carrier.out[hh, yy, ss] for hh in 1:mg.parameters.nh)
                end
            end
        elseif a.carrier  isa Hydrogen
            hydrogen[yy,ss] = sum(max(0., power_balance(hh, yy, ss, mg, typeof(Hydrogen()))) for hh in 1:mg.parameters.nh) / sum(a.carrier.out[hh, yy, ss] for hh in 1:mg.parameters.nh)
            for aa in mg.grids
                if aa.carrier isa Hydrogen
                    hydrogen[yy,ss] = hydrogen[yy,ss] - sum(aa.carrier.in[hh, yy, ss] for hh in 1:mg.parameters.nh) / sum(a.carrier.out[hh, yy, ss] for hh in 1:mg.parameters.nh)
                end
            end
        end
    end

    return LPSP(elec, heat, hydrogen)
end

mutable struct Metrics{T}
    baseline::T
    npv::NPV{T}
    eac::EAC{T}
    renewable_share::T
    lpsp::LPSP{T}
end

# Compute indicators
Metrics(mg::Microgrid, designer::AbstractDesigner) = Metrics(1:mg.parameters.ns, mg, designer)
# Compute indicators for a given scenario s
function Metrics(s::Union{Int64, UnitRange{Int64}}, mg::Microgrid, designer::AbstractDesigner)
    # Baseline cost
    baseline = dropdims(baseline_cost(mg), dims = 1)
    # NPV
    npv = NPV(s, mg, designer)
    # EAC
    eac = EAC(1, s, mg, designer)
    # Share of renewables
    share = renewable_share(s, mg)
    # LPSP
    lpsp = LPSP(s, mg)

    return Metrics(baseline, npv, eac, share, lpsp)
end
