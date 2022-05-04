module CustomSimulation

using NBodySimulator
using Unitful

import NBodySimulator.NBodySimulation
import NBodySimulator.NBodySystem
import NBodySimulator.BoundaryConditions
import NBodySimulator.PotentialNBodySystem
import NBodySimulator.run_simulation
import NBodySimulator.calculate_simulation

export COCOSimulation


struct COCOSimulation{sType <: NBodySystem,
		      bType <: BoundaryConditions,
		      tType <: Unitful.Time}

  system::sType
  tspan::Tuple{tType, tType}
  boundary_conditions::bType
end


function NBodySimulation(system::PotentialNBodySystem, tspan::Tuple{tType,tType}) where {tType <: Unitful.Time}
  COCOSimulation(system, tspan, InfiniteBox())
end

function run_simulation(s::COCOSimulation, args...; kwargs...)
  calculate_simulation(s, args...; kwargs...)
end

function

end
