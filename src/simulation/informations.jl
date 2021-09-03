#=
    This file includes all the funtions needed to update the operation
    and investment informations
 =#


function update_operation_informations!(h::Int64, y::Int64, s::Int64, mg::Microgrid, ω::Scenarios)
    # Demands
    for a in mg.demands
        if a.carrier isa Electricity
            a.timestamp[h,y,s] = ω.ld_E.t[h,y,s]
            a.carrier.out[h,y,s] = ω.ld_E.power[h,y,s]
        elseif d.carrier isa Heat
            a.timestamp[h,y,s] = ω.ld_H.t[h,y,s]
            a.carrier.out[h,y,s] = ω.ld_H.power[h,y,s]
        elseif d.carrier isa Hydrogen
            a.timestamp[h,y,s] = ω.ld_H2.t[h,y,s]
            a.carrier.out[h,y,s] = ω.ld_H2.power[h,y,s]
        end
    end
    # Generations
    for a in mg.generations
        if a.carrier isa Electricity
            a.timestamp[h,y,s] = ω.pv.t[h,y,s]
            a.carrier.in[h,y,s] = a.powerMax[y,s] * ω.pv.power[h,y,s]
        elseif a.carrier isa Heat
            println("Not yet implemented")
        elseif a.carrier isa Hydrogen
            println("Not yet implemented")
        end
    end
    # Grids - We assume the price of electricity is known over the year
    for a in mg.grids
        if a.carrier isa Electricity && h == 1
            a.cost_in[1:end,y,s] = ω.grid.cost_in[1:end,y,s]
            a.cost_out[1:end,y,s] = ω.grid.cost_out[1:end,y,s]
        elseif a.carrier isa Heat && h == 1
            println("Not yet implemented")
        elseif a.carrier isa Hydrogen && h == 1
            println("Not yet implemented")
        end
    end
end
function update_investment_informations!(y::Int64, s::Int64, mg::Microgrid, ω::Scenarios)
    # Generation
    for a in mg.generations
        a.cost[y,s] = ω.pv.cost[y,s]
    end
    # Liion
    for a in mg.storages
        a.cost[y,s] = ω.liion.cost[y,s]
    end
end
