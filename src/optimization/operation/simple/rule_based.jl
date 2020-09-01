#=
    Rule based controller
=#

#                                   Policy definition
#_______________________________________________________________________________
function rb_operation_policy(h::Int64, y::Int64, s::Int64, ld::Load, pv::Source, liion::Liion)
    # Liion
    u_liion = ld.power_E[h,y,s] - pv.power_E[h,y,s]
    return u_liion
end

#                                   Offline functions
#_______________________________________________________________________________
function initialize_controller(ld::Load, pv::Source, liion::Liion,
    controller::RuleBasedController, grid::Grid, ω_optim::Scenarios, parameters::NamedTuple)
     # Parameters
     nh = size(ld.power_E,1) # number of simulation hours in one year
     ny = size(ld.power_E,2) # number of simulation years
     ns = size(ld.power_E,3) # number of scenarios

     # Initialize controller policy
     controller.π = rb_operation_policy

     # Initialize decisions variables
     controller.u = (
     u_liion = convert(SharedArray,zeros(nh,ny,ns)),
     )
end

#                                   Online functions
#_______________________________________________________________________________
function compute_operation_decisions(h::Int64, y::Int64, s::Int64, ld::Load, pv::Source,
     liion::Liion, grid::Grid, controller::RuleBasedController, ω_optim::Scenarios, parameters::NamedTuple)
    controller.u.u_liion[h,y,s] = controller.π(h, y, s, ld, pv, liion)
end
