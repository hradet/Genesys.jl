#=
    This file includes all the abstract struct needed to run the simulation
=#

abstract type AbstractController end
abstract type AbstractDesigner end
#TODO a definir ailleurs ?
abstract type AbstractOneStageDesigner <: AbstractDesigner end
abstract type AbstractOneStageStochasticDesigner <: AbstractDesigner end
abstract type AbstractMultiStageDesigner <: AbstractDesigner end
abstract type AbstractMultiStageStochasticDesigner <: AbstractDesigner end
abstract type AbstractScenarios end
