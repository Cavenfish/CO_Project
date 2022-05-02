module MvH_Potential

using NBodySimulator
using StaticArrays
import NBodySimulator.get_accelerating_function

export MorseParameters
export DispersionParameters

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

struct DispersionParameters <: PotentialParameters
  Ccc::Float64
  Coo::Float64
  Cco::Float64
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

end #module
