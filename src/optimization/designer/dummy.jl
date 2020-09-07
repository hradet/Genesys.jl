#=
    Dummy designer
=#

mutable struct DummyDesigner <: AbstractDesigner
    u::NamedTuple
    DummyDesigner() = new()
end

#### Offline functions ####
# Simple
function initialize_designer(ld::Genesys.Load, pv::Genesys.Source,
    liion::Genesys.Liion, designer::DummyDesigner,
     grid::Genesys.Grid, ω_optim::Genesys.Scenarios, parameters::NamedTuple)
     # Parameters
     nh = size(ld.power_E,1) # number of simulation hours in one year
     ny = size(ld.power_E,2) # number of simulation years
     ns = size(ld.power_E,3) # number of scenarios

     # Initialize controller and designer policies
     # The policy must be initialize at this place...

     # Initialize decisions variables
     designer.u = (
     u_liion = zeros(ny,ns),
     u_pv = zeros(ny,ns),
     )
end

#### Online functions ####
# Simple
function compute_investment_decisions(y::Int64, s::Int64, ld::Load, pv::Source,
    liion::Liion, grid::Grid, designer::DummyDesigner, ω_optim::Scenarios, parameters::NamedTuple)
    ϵ = 0.1
    if liion.soh[end,y,s] < ϵ
        designer.u.u_liion[y,s] = liion.Erated[y,s]
    end
end
# Multi-energy
function compute_investment_decisions(y::Int64, s::Int64, ld::Load, pv::Source,
    liion::Liion, h2tank::H2Tank, elyz::Electrolyzer, fc::FuelCell, tes::ThermalSto,
    heater::Heater, designer::DummyDesigner, ω_optim::Scenarios, parameters::NamedTuple)
    ϵ = 0.1

    # Liion
    if liion.soh[end,y,s] < ϵ
        designer.u.u_liion[y,s] = liion.Erated[y,s]
    end

    # Electrolyzer
    if elyz.soh[end,y,s] < ϵ
        designer.u.u_elyz[y,s] = elyz.powerMax[y,s]
    end

    # FuelCell
    if fc.soh[end,y,s] < ϵ
        designer.u.u_fc[y,s] = fc.powerMax[y,s]
    end
end
