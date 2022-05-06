module CustomSimulation

using Unitful
using UnitfulAtomic

export Vacuo
export MDSimulation

abstract type Simulation
end

abstract type SimulationCell
end

struct InfiniteCell{} <: SimulationCell
  bounds::SVector{6, Real}
end

Vacuo() = InfiniteCell(SVector{-Inf, Inf,-Inf, Inf,-Inf, Inf})

struct MDSimulation{sType <: System,
		    bType <: SimulationCell,
		    tType <: Unitful.Time} <: Simulation

  system::sType
  tspan::Tuple{tType, tType}
  boundary_conditions::bType
end

function MDSimulation(sys::System, tspan::Tuple{tType, tType}, 
		      cell::SimulationCell) where {tType <: Unitful.Time}
  MDSimulation(sys, tspan, cell)
end

end
