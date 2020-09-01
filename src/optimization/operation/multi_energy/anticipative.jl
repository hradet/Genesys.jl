#=
    Anticipative controller
=#

#                                   Model definition
#_______________________________________________________________________________

#                                   Offline functions
#_______________________________________________________________________________

#                                   Online functions
#_______________________________________________________________________________
function compute_operation_decisions(h::Int64, y::Int64, s::Int64, ld::Load, pv::Source,
     liion::Liion, h2tank::H2Tank, elyz::Electrolyzer, fc::FuelCell, tes::ThermalSto,
     heater::Heater, controller::AnticipativeController, ω_optim::Scenarios, parameters::NamedTuple)
    return nothing
end
