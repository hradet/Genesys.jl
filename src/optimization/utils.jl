function preallocate!(controller::AbstractController, nh::Int64, ny::Int64, ns::Int64)
    controller.u = (liion = zeros(nh,ny,ns),
                    elyz = zeros(nh,ny,ns),
                    fc = zeros(nh,ny,ns),
                    tes = zeros(nh,ny,ns),
                    heater = zeros(nh,ny,ns))
end

function preallocate!(designer::AbstractDesigner, ny::Int64, ns::Int64)
    designer.u = (pv = zeros(ny,ns),
                  liion = zeros(ny,ns),
                  tes = zeros(ny,ns),
                  h2tank = zeros(ny,ns),
                  elyz = zeros(ny,ns),
                  fc = zeros(ny,ns))
end
