#=
    Rule based designer with rule based controller
=#

#                                   Policy definition
#_______________________________________________________________________________
function rb_investment_policy(y::Int64, s::Int64, ld::Load, pv::Source, liion::Liion,
      h2tank::H2Tank, elyz::Electrolyzer, fc::FuelCell, tes::ThermalSto, heater::Heater)
    ϵ = 0.01
    # Replacement if soh < ϵ
    r_liion = liion.soh[end,y,s] < ϵ || y==1
    r_elyz = elyz.soh[end,y,s] < ϵ || y==1
    r_fc = fc.soh[end,y,s]< ϵ || y==1
    # First year investment
    δ_inv = y==1

    # Outputs
    u_pv = 50 * δ_inv
    u_liion = 200. * r_liion
    u_tank = 20000. * δ_inv
    u_elyz = 10. * r_elyz
    u_fc = 5. * r_fc
    u_tes = 200. * δ_inv

    return u_pv, u_liion, u_tank, u_elyz, u_fc, u_tes
end

#                                   Offline functions
#_______________________________________________________________________________
function initialize_designer(ld::Load, pv::Source, liion::Liion, h2tank::H2Tank,
      elyz::Electrolyzer, fc::FuelCell, tes::ThermalSto, heater::Heater,
      designer::RuleBasedDesigner, grid::Grid, ω_optim::Scenarios, parameters::NamedTuple)
      # Parameters
      nh = size(ld.power_E,1) # number of simulation hours in one year
      ny = size(ld.power_E,2) # number of simulation years
      ns = size(ld.power_E,3) # number of scenarios

      # Initialize designer policy
      designer.π = rb_investment_policy

      # Initialize decisions variables
      designer.u = (
      u_liion = zeros(ny,ns),
      u_pv = zeros(ny,ns),
      u_tank = zeros(ny,ns),
      u_elyz = zeros(ny,ns),
      u_fc = zeros(ny,ns),
      u_tes = zeros(ny,ns),
      )
end

#                                   Online functions
#_______________________________________________________________________________
function compute_investment_decisions(y::Int64, s::Int64, ld::Load, pv::Source,
    liion::Liion, h2tank::H2Tank, elyz::Electrolyzer, fc::FuelCell, tes::ThermalSto,
    heater::Heater, designer::RuleBasedDesigner, ω_optim::Scenarios, parameters::NamedTuple)
    designer.u.u_pv[y,s], designer.u.u_liion[y,s], designer.u.u_tank[y,s], designer.u.u_elyz[y,s], designer.u.u_fc[y,s], designer.u.u_tes[y,s] = designer.π(y, s, ld, pv, liion, h2tank, elyz, fc, tes, heater)
end
