# Reshape from time block to simulation for operation variables
function reshape_tb_to_simu_operation(u_opt, n_bytb; nh=8760)
    # Parameters
    u_simu = []
    ntb = size(u_opt,2) # number of time blocks
    nhh = size(u_opt,1) # number of hours
    # Append u_simu with repeated values of u_opt
    for tb in 1:ntb
        append!(u_simu, repeat(u_opt[:,tb], ceil(Int,nh/nhh) * n_bytb[tb],1))
    end
    # Formatting
    u_simu = reshape(u_simu, nh, :)
    return u_simu
end
# Reshape from time block to simulation for investment variables
function reshape_tb_to_simu_investment(u_opt, σ_tb)
    # Pre allocate
    u_simu = zeros(σ_tb[end])
    # Replace by u_opt at the given year
    u_simu[σ_tb] = u_opt
    return u_simu
end
