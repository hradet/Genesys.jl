#=
    This file includes all the funtions needed to compute the operation
    and investment dynamics
 =#

function compute_operation_dynamics!(h::Int64, y::Int64, s::Int64, mg::Microgrid, controller::AbstractController)
    # Converters
    for (k, a) in enumerate(mg.converters)
        # Attention pas si simple car plusieurs carrier pour chaque converter...
        a.carrier.in[h,y,s], a.carrier.out[h,y,s], a.soh[h+1,y,s] =
        compute_operation_dynamics(a, (powerMax = a.powerMax[y,s], soh = a.soh[h,y,s]), controller.decisions.converters[k][h,y,s], mg.parameters.Δh)
    end
    # Storage
    for (k, a) in enumerate(mg.storages)
        a.soc[h+1,y,s], a.soh[h+1,y,s], a.carrier.in[h,y,s], a.carrier.out[h,y,s] =
        compute_operation_dynamics(a, (Erated = a.Erated[y,s], soc = a.soc[h,y,s], soh = a.soh[h,y,s]), controller.decisions.storages[k][h,y,s], mg.parameters.Δh)
    end
end
function compute_investment_dynamics!(y::Int64, s::Int64, mg::Microgrid, designer::AbstractDesigner)
     # Generations
     for (k, a) in enumerate(mg.generations)
         a.powerMax[y+1,s] = compute_investment_dynamics(a, (powerMax = a.powerMax[y,s],), designer.decisions.generations[k][y,s])
     end
     # Storages
     for (k, a) in enumerate(mg.storages)
         a.Erated[y+1,s], a.soc[1,y+1,s], a.soh[1,y+1,s] =
         compute_investment_dynamics(a, (Erated = a.Erated[y,s], soc = a.soc[end,y,s], soh = a.soh[end,y,s]), designer.decisions.storages[k][y,s])
     end
     # Converters
     for (k, a) in enumerate(mg.converters)
         a.powerMax[y+1,s], a.soh[1,y+1,s] =
         compute_investment_dynamics(a, (powerMax = a.powerMax[y,s], soh = a.soh[end,y,s]), designer.decisions.converters[k][y,s])
     end
end
