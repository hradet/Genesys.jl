#=
    Abstract designer
=#

abstract type AbstractDesigner end

# Preallocate abstract designer
function preallocate!(mg::Microgrid, designer::AbstractDesigner)
    designer.decisions = (generations = [zeros(mg.parameters.ny, mg.parameters.ns) for a in mg.generations],
                          storages = [zeros(mg.parameters.ny, mg.parameters.ns) for a in mg.storages],
                          converters = [zeros(mg.parameters.ny, mg.parameters.ns) for a in mg.converters])
end
