include("./AutoMD.jl")
include("MvH_Potential.jl")

using .AutoMD
using .MvH_Potential
using NBodySimulator
using StaticArrays

l    = 10.0
bdys = read_ase_xyz("10co.xyz")
#pbc  = PeriodicBoundaryConditions(l)
pbc  = InfiniteBox()
#sys  = ChargedParticles(bdys, 1.0)

p = MorseParameters(11.23, 2.3281*1.1282, 1.1282)
sys = PotentialNBodySystem(bdys, Dict(:custom_potential_params => p))

sim  = NBodySimulation(sys, (0.0, 1.0), pbc)
sr   = run_simulation(sim, VelocityVerlet(), dt=0.001)
#NBodySimulator.save_to_pdb(sr, "test.pdb")

write_xyz_traj("test.xyz", sr)
