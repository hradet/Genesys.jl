#=
    Rule based designer with rule based controller
=#

#                                   Policy definition
#_______________________________________________________________________________
function rb_investment_policy(y::Int64, s::Int64, ld::Load, pv::Source, liion::Liion)
    ϵ = 0.01
    # Replacement if soh < ϵ or first year
    δ_inv = liion.soh[end,y,s] < ϵ || y==1
    # Daily power net
    daily_power_net=copy(reshape(ld.power_E[1:end,y,s], 24, :))
    # We only keep positive values
    daily_power_net[daily_power_net .< 0.] .= 0.
    # Mean daily energy consumption
    mean_daily = mean(sum(daily_power_net,dims =1))
    max_daily = maximum(daily_power_net)
    # Outputs
    u_liion = mean_daily / 2 * δ_inv
    u_pv = 2 * max_daily * δ_inv

    return u_pv, u_liion
end

#                                   Offline functions
#_______________________________________________________________________________
function initialize_designer(ld::Load, pv::Source, liion::Liion,
     designer::RuleBasedDesigner, grid::Grid, ω_optim::Scenarios, parameters::NamedTuple)
     # Parameters
     nh = size(ld.power_E,1) # number of simulation hours in one year
     ny = size(ld.power_E,2) # number of simulation years
     ns = size(ld.power_E,3) # number of scenarios

     # Initialize designer policy
     designer.π = rb_investment_policy

     # Initialize decisions variables
     designer.u = (
     u_liion = convert(SharedArray,zeros(ny,ns)),
     u_pv = convert(SharedArray,zeros(ny,ns)),
     )
end

#                                   Online functions
#_______________________________________________________________________________
function compute_investment_decisions(y::Int64, s::Int64, ld::Load, pv::Source,
     liion::Liion, grid::Grid, controller::AbstractController,
     designer::RuleBasedDesigner, ω_optim::Scenarios, parameters::NamedTuple)
    designer.u.u_pv[y,s], designer.u.u_liion[y,s] = designer.π(y, s, ld, pv, liion)
end
