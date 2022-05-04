module ExchangePotential

include("CustomBodies.jl")

using .CustomBodies
using NBodySimulator
using StaticArrays
using Unitful
using UnitfulAtomic
import NBodySimulator.get_accelerating_function


export ExchangeParameters

struct ExchangeParameters <: PotentialParameters
  Acc::typeof(1.0u"eV")
  Aoo::typeof(1.0u"eV")
  Aco::typeof(1.0u"eV")
  Bcc::typeof(1.0u"Ang^-1")
  Boo::typeof(1.0u"Ang^-1")
  Bco::typeof(1.0u"Ang^-1")
end

function get_accelerating_function(p::ExchangeParameters, sim::NBodySimulation)
  #Note, dv = acceleration
  #      u  = positions
  #      v  = velocity
  #      t  = time
  #      i  = particle number

  function V_exch(diff, Aij, Bij)
    r = √(diff'diff)
    E = Aij * exp(-Bij * r)
    F = Bij * E * diff / r
    return F
  end

  bdy = sim.system.bodies
  m   = get_masses(sim.system)
  N   = length(bdy)

  (dv, u, v, t, i) -> begin

    F  = @SVector [0.0, 0.0, 0.0]
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

      if     (bdy[i].s == bdy[j].s == 'C')
        F += V_exch(rj - ri, p.Acc, p.Bcc)
      elseif (bdy[i].s == bdy[j].s == 'O')
        F += V_exch(rj - ri, p.Aoo, p.Boo)
      else
        F += V_exch(rj - ri, p.Aco, p.Bco)
      end
    end

    dv .= F / m[i]
  end # Exchange acceleration
end

end
