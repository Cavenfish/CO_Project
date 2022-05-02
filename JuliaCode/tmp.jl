include("./AutoMD.jl")
include("MvH_Potential.jl")

using .AutoMD
using .MvH_Potential
using NBodySimulator
using StaticArrays

l    = 10.0
bdys = read_ase_xyz("200co5_NVE.xyz")
#pbc  = PeriodicBoundaryConditions(l)
pbc  = InfiniteBox()
#sys  = ChargedParticles(bdys, 1.0)

M = MorseParameters(11.23, 2.3281*1.1282, 1.1282)
d = DispersionParameters(-33.37, -10.52, -15.16)

pot = Dict(:custom_potential_params => M,
           :custom_potential_params => d)

sys = PotentialNBodySystem(bdys, pot)

τ    = 1e-5
sim  = NBodySimulation(sys, (0.0, 1000τ), pbc)
sr   = run_simulation(sim, VelocityVerlet(), dt=τ)
#NBodySimulator.save_to_pdb(sr, "test.pdb")

write_xyz_traj("test.xyz", sr, 1)
