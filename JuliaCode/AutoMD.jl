module AutoMD

include("CustomBodies.jl")

using .CustomBodies
using StaticArrays
using NBodySimulator
using Unitful
using UnitfulAtomic

export read_ase_xyz
export write_xyz_traj

function read_ase_xyz(xyz)
  sys = readlines(xyz)
  N   = length(sys) - 2
  hed = sys[2]
  sys = split.(sys[3:end], " ")
  sys = deleteat!.(sys, findall.(e -> e == "", sys))
  amu = Dict("C" => 12.011u"u", "O" => 15.999u"u")
  Q   = Dict("C" => -1.786u"e_au", "O" => -2.332u"e_au")
  set = Nucleus[]

  for i in range(1,N)
    props = parse.(Float64, sys[i][2:end])
    pos   = SVector{3}(props[1:3])u"Ang"
    s     = sys[i][1]
    q     = Q[s]

    if occursin("masses", hed)
      mas = props[4]
      vel = SVector{3}(props[5:7]u"(eV*u)^0.5" ./ mas)
    else
      mas = amu[s]
      vel = SVector{3}(props[4:6]u"(eV*u)^0.5" ./ mas)
    end #if-else

    particle = Nucleus(pos, vel, mas, q, s)
    push!(set, particle)
  end #for loop

  return set
end #read_ase_xyz


function write_xyz_traj(fileName::String, sr, dt)
  bodies = sr.simulation.system.bodies
  n      = length(bodies)
  i      = 0
  f      = open(fileName, "w")

  for t in sr.solution.t[1:dt:end]
    cc = get_position(sr, t)
    i += 1

    println(f, lpad(n, 9))
    println(f, "i=$i, time=$t")

    for j ∈ 1:n

      s = bodies[j].s
      x = cc[1,j]
      y = cc[2,j]
      z = cc[3,j]

      println(f, "$s   $x   $y   $z")
    end # molecule for loop

  end # time for loop
  close(f)
end #write_xyz_traj

end #module
