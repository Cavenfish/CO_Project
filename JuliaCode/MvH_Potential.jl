module MvH_Potential

using NBodySimulator
using StaticArrays
import NBodySimulator.get_accelerating_function

export MorseParameters

struct MorseParameters <: PotentialParameters
  ϵ::Float64
  ρ0::Float64
  r0::Float64
  # ϵ  = 11.23
  # ρ0 = 2.3281*1.1282
  # r0 = 1.1282
end

function get_accelerating_function(p::MorseParameters, sim::NBodySimulation)
  #Note, dv = acceleration
  #      u  = positions
  #      v  = velocity
  #      t  = time
  #      i  = particle number

  bodies = sim.system.bodies
  m = get_masses(sim.system)
  N      = length(bodies)

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
    F   += preF * expf * (expf - 1) * diff / r

    dv  .= F / m[i]
  end #Morse acceleration
end

end #module
