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
  r0::Float64
  FQ::Bool
  FX::Bool
end

function get_accelerating_function(p::CoulombParameters, sim::NBodySimulation)
  #Note, dv = acceleration
  #      u  = positions
  #      v  = velocity
  #      t  = time
  #      i  = particle number

  function FQ(ri, rj, qi, qj, α, v)
    d = rj - ri
    r = √(d'd)
    E = qi * qj / r
    F = α * E * v
    return F
  end

  function FQX(ri, rj, qi, qj, c, qx, v)
    d = rj - ri
    r = √(d'd)
    E = qi * qj / r
    F = c * E / qx * v
    return F
  end

  function F_coul(diff, Qij)
    r = √(diff'diff)
    F = Qij * diff / r^3
    return F
  end

  function dist(v1, v2)
    diff = v1 .- v2
    r    = √(diff'diff)
    return r
  end

  bdy = sim.system.bodies
  m   = get_masses(sim.system)
  N   = length(bdy)
  W   = Dict("C" => 0.57135, "O" => 0.42865)
  Q   = Dict("C" => p.αc,    "O" => p.αo)

  (dv, u, v, t, i) -> begin
    ib  = bdy[i].b
    F   = @SVector [0.0, 0.0, 0.0]
    ri  = @SVector [u[1,i ], u[2,i ], u[3,i ]]
    rib = @SVector [u[1,ib], u[2,ib], u[3,ib]]
    di  = dist(ri, rib)
    qi  = bdy[i].q * exp(-Q[bdy[i].s] * (di - p.r0))

    for j in 1:N

      if (i ==  j) || (j == bdy[i].b)
        continue
      end

      jb  = bdy[j].b
      rj  = @SVector [u[1,j ], u[2,j ], u[3,j ]]
      rjb = @SVector [u[1,jb], u[2,jb], u[3,jb]]
      dj  = dist(rj, rjb)
      qj  = bdy[j].q * exp(-Q[bdy[j].s] * (dj - p.r0))

      F -= F_coul(rj .- ri, qi * qj)

      if p.FX
        qib = bdy[ib].q * exp(-Q[bdy[ib].s] * (di - p.r0))
        qjb = bdy[jb].q * exp(-Q[bdy[jb].s] * (dj - p.r0))
        rXi = (W[bdy[i].s] * ri) + (W[bdy[ib].s] * rib)
        rXj = (W[bdy[j].s] * rj) + (W[bdy[jb].s] * rjb)
        qXi = - (qi + qib)
        qXj = - (qj + qjb)

        F -= F_coul(rXj -  ri, qXj *  qi) * 0.5
        F -= F_coul( rj - rXi,  qj * qXi) * W[bdy[i].s]
        F -= F_coul(rXj - rXi, qXj * qXi) * W[bdy[i].s] * 0.5
      end

      if p.FQ
        #FQ(ri, rj, qi, qj, α, v)
        #FQX(ri, rj, qi, qj, c, qx, v)
        R = rib - ri
        v = R / di
        c = - (Q[bdy[i].s] * qi + Q[bdy[ib].s] * qib)

        F -= FQ(ri, rj, qi, qj, Q[bdy[i].s], v)
        F -= FQ(ri, rXj, qi, qXj, Q[bdy[i].s], v) * 0.5
        F -= FQ(ri, rjb, qi, qjb, Q[bdy[i].s], v)
        F -= FQ(rib, rj, qib, qj, Q[bdy[ib].s], v)
        F -= FQ(rib, rXj, qib, qXj, Q[bdy[ib].s], v) * 0.5
        F -= FQ(rib, rjb, qib, qjb, Q[bdy[ib].s], v)

        F -= FQX(rXi, rj, qXi, qj, c, qXi, v)
        F -= FQX(rXi, rXj, qXi, qXj, c, qXi, v) * 0.5
        F -= FQX(rXi, rjb, qXi, qjb, c, qXi, v)
      end
    end

    dv  .+= F / m[i]
  end # Coulomb acceleration
end

end
