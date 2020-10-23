#=
    This file includes all the funtions needed to compute the techno-economic
    indicators
 =#


 mutable struct Costs
      capex::Union{Array{Float64,1}, Array{Float64,2}}
      opex::Union{Array{Float64,1}, Array{Float64,2}}
      cf::Union{Array{Float64,1}, Array{Float64,2}}
      cumulative_npv::Union{Array{Float64,1}, Array{Float64,2}}
      npv::Union{Float64, Array{Float64,1}}
 end

 # Compute costs
function Costs(des::DistributedEnergySystem, designer::AbstractDesigner)

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
# Compute costs for a given scenario s
function Costs(s::Int64, des::DistributedEnergySystem, designer::AbstractDesigner)

     # Discount factor
     γ = 1. ./ (1. + des.parameters.τ) .^ range(1, length = des.parameters.ny, step = des.parameters.Δy)

    # Discounted capex
    capex = γ .* compute_capex(s, des, designer)

    # Discounted opex
    opex = γ .* compute_opex(s, des)

    # Discounted cash flow
    cf = - capex .+ opex

    # Discounted NPV each year
    cumulative_npv = cumsum(cf, dims = 1)

    # Discounted NPV
    npv = sum(cf, dims=1)

    return Costs(capex, opex, cf, cumulative_npv, npv)
end

# OPEX relative to the baseline cost
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
# OPEX relative to the baseline cost for a given scenario s
function compute_opex(s::Int64, des::DistributedEnergySystem)
    # Total load
    ld_tot = zeros(des.parameters.nh, des.parameters.ny)
    isa(des.ld_E, Load) ? ld_tot .+= des.ld_E.power[:,:,s] : nothing
    isa(des.ld_H, Load) ? ld_tot .+= des.ld_H.power[:,:,s] ./ des.heater.η_E_H : nothing

    # Reference case when all the electricity is purchased from the grid
    ref = max.(0, ld_tot)
    saving = ((ref .- max.(0., des.grid.power_E[:,:,s])) .* des.grid.C_grid_in[:,:,s] - min.(0., des.grid.power_E[:,:,s]) .* des.grid.C_grid_out[:,:,s]) .* des.parameters.Δh

    # opex
    opex = dropdims(sum(saving, dims=1), dims=1)

    return opex
end
# OPEX for eac
function compute_opex_eac(y::Int64, s::Int64, des::DistributedEnergySystem)
    opex =  sum(max.(0., des.grid.power_E[:,y,s]) .* des.grid.C_grid_in[:,y,s] - min.(0., des.grid.power_E[:,y,s]) .* des.grid.C_grid_out[:,y,s]) * des.parameters.Δh
    return opex
end

# CAPEX
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
# CAPEX for a given scenario s
function compute_capex(s::Int64, des::DistributedEnergySystem, designer::AbstractDesigner)
    # Preallocation
    capex = zeros(des.parameters.ny)

    isa(des.pv, Source) ? capex .+= designer.u.pv[:,s] .* des.pv.C_pv[:,s] : nothing
    isa(des.liion, Liion) ? capex .+= designer.u.liion[:,s] .* des.liion.C_liion[:,s] : nothing
    isa(des.tes, ThermalSto) ? capex .+= designer.u.tes[:,s] .* des.tes.C_tes[:,s] : nothing
    isa(des.h2tank, H2Tank) ? capex .+= designer.u.h2tank[:,s] .* des.h2tank.C_tank[:,s] : nothing
    isa(des.elyz, Electrolyzer) ? capex .+= designer.u.elyz[:,s] .* des.elyz.C_elyz[:,s] : nothing
    isa(des.fc, FuelCell) ? capex .+= designer.u.fc[:,s] .* des.fc.C_fc[:,s] : nothing

    return capex
end
# CAPEX for eac
function compute_capex_eac(y::Int64, s::Int64, des::DistributedEnergySystem, designer::AbstractDesigner)
    # Preallocation
    capex = 0.

    if isa(des.pv, Source)
        Γ_pv = (des.parameters.τ * (des.parameters.τ + 1.) ^ des.pv.lifetime) / ((des.parameters.τ + 1.) ^ des.pv.lifetime - 1.)
        capex += Γ_pv * designer.u.pv[y,s] * des.pv.C_pv[y,s]
    end

    if isa(des.liion, Liion)
        Γ_liion = (des.parameters.τ * (des.parameters.τ + 1.) ^ des.liion.lifetime) / ((des.parameters.τ + 1.) ^ des.liion.lifetime - 1.)
        capex += Γ_liion * designer.u.liion[y,s] * des.liion.C_liion[y,s]
    end

    if isa(des.tes, ThermalSto)
        Γ_tes = (des.parameters.τ * (des.parameters.τ + 1.) ^ des.tes.lifetime) / ((des.parameters.τ + 1.) ^ des.tes.lifetime - 1.)
        capex += Γ_tes * designer.u.tes[y,s] * des.tes.C_tes[y,s]
    end

    if isa(des.h2tank, H2Tank)
        Γ_h2tank = (des.parameters.τ * (des.parameters.τ + 1.) ^ des.h2tank.lifetime) / ((des.parameters.τ + 1.) ^ des.h2tank.lifetime - 1.)
        capex += Γ_h2tank * designer.u.h2tank[y,s] * des.h2tank.C_tank[y,s]
    end

    if isa(des.elyz, Electrolyzer)
        Γ_elyz= (des.parameters.τ * (des.parameters.τ + 1.) ^ des.elyz.lifetime) / ((des.parameters.τ + 1.) ^ des.elyz.lifetime - 1.)
        capex += Γ_elyz * designer.u.elyz[y,s] * des.elyz.C_elyz[y,s]
    end

    if isa(des.fc, FuelCell)
        Γ_fc = (des.parameters.τ * (des.parameters.τ + 1.) ^ des.fc.lifetime) / ((des.parameters.τ + 1.) ^ des.fc.lifetime - 1.)
        capex += Γ_fc * designer.u.fc[y,s] * des.fc.C_fc[y,s]
    end

    return capex
end

# Share of renewables
function compute_share(des::DistributedEnergySystem)
    # Total demand
    ld_tot = zeros(des.parameters.ny, des.parameters.ns)
    isa(des.ld_E, Load) ? ld_tot .+= dropdims(sum(des.ld_E.power, dims=1), dims=1) : nothing
    isa(des.ld_H, Load) ? ld_tot .+= dropdims(sum(des.ld_H.power ./ des.heater.η_E_H, dims=1), dims=1) : nothing

    # Share of renew.
    τ_share = 1. .- dropdims(sum(max.(0., des.grid.power_E), dims=1), dims=1) ./ ld_tot

    return τ_share
end
# Share of renewables for a given scenario s
function compute_share(s::Int64, des::DistributedEnergySystem)
    # Total demand
    ld_tot = zeros(des.parameters.ny, des.parameters.ns)
    isa(des.ld_E, Load) ? ld_tot .+= dropdims(sum(des.ld_E.power[:,:,s], dims=1), dims=1) : nothing
    isa(des.ld_H, Load) ? ld_tot .+= dropdims(sum(des.ld_H.power[:,:,s] ./ des.heater.η_E_H, dims=1), dims=1) : nothing

    # Share of renew.
    τ_share = 1. .- dropdims(sum(max.(0., des.grid.power_E[:,:,s]), dims=1), dims=1) ./ ld_tot

    return τ_share
end
# Share of renewables for a given year y of a givn scenario s
function compute_share(y::Int64, s::Int64, des::DistributedEnergySystem)
    # Total demand
    ld_tot = 0.
    isa(des.ld_E, Load) ? ld_tot += sum(des.ld_E.power[:,y,s]) : nothing
    isa(des.ld_H, Load) ? ld_tot += sum(des.ld_H.power[:,y,s]) ./ des.heater.η_E_H : nothing

    # Share of renew.
    τ_share = 1. .- sum(max.(0., des.grid.power_E[:,y,s])) ./ ld_tot

    return τ_share
end

mutable struct LPSP
    lpsp_E::Union{Nothing, Array{Float64,1}, Array{Float64,2}}
    lpsp_H::Union{Nothing, Array{Float64,1}, Array{Float64,2}}
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

# Equivalent annual cost
function compute_eac(y::Int64, s::Int64, des::DistributedEnergySystem, designer::AbstractDesigner)
    # CAPEX
    capex = compute_capex_eac(y, s, des, designer)

    # OPEX
    opex = compute_opex_eac(y, s, des)

    return capex + opex
end

mutable struct Metrics
    costs::Costs
    τ_share::Union{Array{Float64,1}, Array{Float64,2}}
    lpsp::LPSP
end

# Compute indicators
function Metrics(des::DistributedEnergySystem, designer::AbstractDesigner)
    ### Econmics
    costs = Costs(des, designer)

    ### Share of renewables
    τ_share = compute_share(des)

    ### LPSP
    lpsp = LPSP(des)

    return Metrics(costs, τ_share, lpsp)
end
# Compute indicators for a given scenario s
function Metrics(s::Int64, des::DistributedEnergySystem, designer::AbstractDesigner)
    ### Econmics
    costs = Costs(s, des, designer)

    ### Share of renewables
    τ_share = compute_share(s, des)

    ### LPSP
    lpsp = LPSP(s, des)

    return Metrics(costs, τ_share, lpsp)
end
