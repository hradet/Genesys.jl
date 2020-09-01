#=
    This file includes all the funtions needed to compute the operation
    and investment dynamics
 =#

# Compute operation dynamic
function compute_operation_dynamics(h::Int64, y::Int64, s::Int64, liion::Liion, controller::AbstractController, parameters::NamedTuple)
    liion.soc[h+1,y,s], liion.soh[h+1,y,s], liion.power_E[h,y,s] = compute_operation_dynamics(liion, (Erated = liion.Erated[y,s], soc = liion.soc[h,y,s], soh = liion.soh[h,y,s]), controller.u.u_liion[h,y,s], parameters.Δh)
end

# Compute investment dynamic
function compute_investment_dynamics(y::Int64, s::Int64, pv::Source, liion::Liion, designer::AbstractDesigner)
     # Converters
     pv.powerMax[y+1,s] = compute_investment_dynamics(pv, (powerMax = pv.powerMax[y,s],), designer.u.u_pv[y,s])
     # Storage
     liion.Erated[y+1,s], liion.soc[1,y+1,s], liion.soh[1,y+1,s] = compute_investment_dynamics(liion, (Erated = liion.Erated[y,s], soc = liion.soc[end,y,s], soh = liion.soh[end,y,s]), designer.u.u_liion[y,s])
end
