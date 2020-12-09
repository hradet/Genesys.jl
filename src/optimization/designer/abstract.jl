#=
    Abstract designer
=#

abstract type AbstractDesigner end

# Preallocate abstract designer
function preallocate!(designer::AbstractDesigner, ny::Int64, ns::Int64)
    designer.u = (pv = zeros(ny,ns),
                  liion = zeros(ny,ns),
                  tes = zeros(ny,ns),
                  h2tank = zeros(ny,ns),
                  elyz = zeros(ny,ns),
                  fc = zeros(ny,ns))
end
