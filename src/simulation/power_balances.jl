#=
    This file includes all the functions needed to check the power balance
    constraints
=#

#TODO rajouter fonction add_to_powerbalance!(power_balance, power)

function compute_power_balances!(h::Int64, y::Int64, s::Int64, mg::Microgrid)
    # Hydrogen
    checking!(h, y, s, mg, typeof(Hydrogen()))
    # Heat
    checking!(h, y, s, mg, typeof(Heat()))
    # Electricity
    checking!(h, y, s, mg, typeof(Electricity()))
end
function power_balance(h::Union{Int64, UnitRange{Int64}}, y::Union{Int64, UnitRange{Int64}}, s::Union{Int64, UnitRange{Int64}}, mg::Microgrid, type::DataType)
    # Parameters
    ϵ = 0.01 # 0.01 kW tolerance
    # Energy balance
    balance = 0.
    # Demands
    for a in mg.demands
        !isempty(mg.demands) && a.carrier isa type ? balance = balance .+ a.carrier.out[h,y,s] : nothing
    end
    # Generations
    for a in mg.generations
        !isempty(mg.generations) && a.carrier isa type ? balance = balance .- a.carrier.in[h,y,s] : nothing
    end
    # Storages
    for a in mg.storages
        !isempty(mg.storages) && a.carrier isa type ? balance = balance .- a.carrier.in[h,y,s] .- a.carrier.out[h,y,s] : nothing
    end
    # Converters
    for a in mg.converters
        if !isempty(mg.converters)
            for c in a.carrier
                c isa type ? balance = balance .- c.in[h,y,s] .- c.out[h,y,s] : nothing
            end
        end
    end
    return balance
end
function checking!(h::Int64, y::Int64, s::Int64, mg::Microgrid, type::DataType)
    # Parameters
    ϵ = 0.01 # 0.01 kW tolerance
    # Energy balance
    balance = power_balance(h, y, s, mg, type)
    # Grids
    for a in mg.grids
        if !isempty(mg.grids)
            if a.carrier isa type
                a.carrier.in[h,y,s] = min(a.powerMax, max(0, balance))
                a.carrier.out[h,y,s] = max(-a.powerMax, min(0, balance))
            end
        else # If the balance equation is not fulfilled, systems are turned to zerp
            if balance > sum(a.carrier.out[h,y,s] for a in mg.demands if a.carrier isa type) + ϵ
                # Storage set to zero
                for a in mg.storages
                    if a.carrier isa type
                        a.carrier.in[h,y,s], a.carrier.out[h,y,s]= 0.
                        a.soc[h+1,y,s] = max(0., a.soc[h,y,s] * (1. - a.η_self))
                    end
                end
                # Converters set to zero
                for a in mg.converters
                    for c in a.carrier
                        c.in[h,y,s], c.out[h,y,s]= 0., 0.
                    end
                end
            else
                println("Shedding energy demand!")
            end
        end
    end
end
