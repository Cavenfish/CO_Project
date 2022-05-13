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
   ϵ::Float64
  ρ0::Float64
  r0::Float64
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

    j    = bdy[i].b
    rj   = @SVector [u[1, j], u[2, j], u[3, j]]

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
