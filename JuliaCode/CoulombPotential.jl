module CoulombPotential

include("CustomBodies.jl")

using .CustomBodies
using NBodySimulator
using StaticArrays
using Unitful
using UnitfulAtomic
import NBodySimulator.get_accelerating_function

export CoulombParameters

struct CoulombParameters <: PotentialParameters
  αc::Float64
  αo::Float64
  r0::typeof(1.0u"Ang")
  FQ::Bool
  FX::Bool
end

function get_accelerating_function(p::CoulombParameters, sim::NBodySimulation)
  #Note, dv = acceleration
  #      u  = positions
  #      v  = velocity
  #      t  = time
  #      i  = particle number

  function V_coul(diff, Qij)
    r = √(diff'diff)
    F = Qij * diff / r^3
    return F
  end

  bdy = sim.system.bodies
  m   = get_masses(sim.system)
  N   = length(bdy)

  (dv, u, v, t, i) -> begin
    F  = @SVector [0.0, 0.0, 0.0]
    FQ = @SVector [0.0, 0.0, 0.0]
    FX = @SVector [0.0, 0.0, 0.0]
    ri = @SVector [u[1,i], u[2,i], u[3,i]]

    #= The current system for avoinding double counting
    and forcing pairwise is very ugly. This needs to
    improve =#

    for j in 1:N
      if (i ==  j)
        continue
      end

      if (i % 2 == 0)
        if (i - 1 == j)
          continue
        end
      else
        if (i + 1 == j)
          continue
        end
      end

      rj = @SVector [u[1,j], u[2,j], u[3,j]]
      F -= V_coul(rj - ri, bdy[i].q * bdy[j].q)
    end

    #=
    if p.FQ
      F += FQ
    end

    if p.FX
      F += FX
    end
    =#

    dv  .= F / m[i]
  end # Coulomb acceleration 
end

end
