module DispersionPotential

include("CustomBodies.jl")

using .CustomBodies
using NBodySimulator
using StaticArrays
using Unitful
using UnitfulAtomic
import NBodySimulator.get_accelerating_function

export DispersionParameters

struct DispersionParameters <: PotentialParameters
  Ccc::typeof(1.0u"eV*Ang^6")
  Coo::typeof(1.0u"eV*Ang^6")
  Cco::typeof(1.0u"eV*Ang^6")
end

function get_accelerating_function(p::DispersionParameters, sim::NBodySimulation)
  #Note, dv = acceleration
  #      u  = positions
  #      v  = velocity
  #      t  = time
  #      i  = particle number

  function V_disp(diff, Cij)
    r2 = diff'diff
    E  = Cij / r2^3
    F  = 6.0 * E * diff / r2
    return F
  end

  bdy = sim.system.bodies
  m   = get_masses(sim.system)
  N   = length(bdy)

  (dv, u, v, t, i) -> begin

    F  = @SVector [0.0, 0.0, 0.0]
    ri = @SVector [u[1, i], u[2, i], u[3, i]]

    #= The current system for avoinding double counting
    and forcing pairwise is very ugly. This needs to
    improve =#

    for j in 1:N
      if (i == j)
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

      rj   = @SVector [u[1, j], u[2, j], u[3, j]]

      if     (bdy[i].s == bdy[j].s == 'C')
        F += V_disp(rj - ri, p.Ccc)
      elseif (bdy[i].s == bdy[j].s == 'O')
        F += V_disp(rj - ri, p.Coo)
      else
        F += V_disp(rj - ri, p.Cco)
      end
    end

    dv  .= F / m[i]
  end #Dispersion acceleration
end



end