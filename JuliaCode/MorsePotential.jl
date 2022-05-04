module MorsePotential

include("CustomBodies.jl")

using .CustomBodies
using NBodySimulator
using StaticArrays
using Unitful
using UnitfulAtomic
import NBodySimulator.get_accelerating_function

export MorseParameters

struct MorseParameters <: PotentialParameters
   ϵ::typeof(1.0u"eV")
  ρ0::Float64
  r0::typeof(1.0u"Ang")
end


function get_accelerating_function(p::MorseParameters, sim::NBodySimulation)
  #Note, dv = acceleration
  #      u  = positions
  #      v  = velocity
  #      t  = time
  #      i  = particle number

  bdy = sim.system.bodies
  m   = get_masses(sim.system)
  N   = length(bdy)

  (dv, u, v, t, i) -> begin

    if ( i % 2 == 0)
      rj   = @SVector [u[1, i-1], u[2, i-1], u[3, i-1]]
    else
      rj   = @SVector [u[1, i+1], u[2, i+1], u[3, i+1]]
    end

    ri   = @SVector [u[1, i], u[2, i], u[3, i]]
    diff = rj .- ri

    F    = @SVector [0.0, 0.0, 0.0]
    preF = 2 * p.ϵ * p.ρ0 / p.r0
    r    = √(diff'diff)
    expf = exp(p.ρ0 * (1.0 - r / p.r0))
    F   -= preF * expf * (expf - 1) * diff / r

    dv  .= F / m[i]
  end #Morse acceleration
end


end
