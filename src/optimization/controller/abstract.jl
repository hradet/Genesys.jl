#=
    Abstract controller
=#

abstract type AbstractController end

# Preallocate abstract controller
function preallocate!(controller::AbstractController, nh::Int64, ny::Int64, ns::Int64)
    controller.u = (liion = zeros(nh,ny,ns),
                    elyz = zeros(nh,ny,ns),
                    fc = zeros(nh,ny,ns),
                    h2tank = zeros(nh,ny,ns),
                    tes = zeros(nh,ny,ns),
                    heater = zeros(nh,ny,ns))
end
