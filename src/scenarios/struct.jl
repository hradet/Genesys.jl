mutable struct ClusteredScenarios
    clusters::Scenarios # clusters
    σ # sequence of clusters = assignments
    nby # number of entities by cluster = counts
end

mutable struct MarkovChain
    states
    transition_matrices
end
