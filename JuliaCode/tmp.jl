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

M = MorseParameters(11.23, 2.3281*1.1282, 1.1282)
d = DispersionParameters(-33.37, -10.52, -15.16)
c = CoulombParameters(-1.786, -2.332, 3.845, 2.132, 1.1282, false, false)
e = ExchangeParameters(361.36, 6370.10, 1516.74, 2.836, 4.253, 3.544)

pot = Dict(:Morse_params      => M,
	   :Dispersion_params => d,
	   :Coulomb_params    => c,
	   :Exchange_params   => e)

sys = PotentialNBodySystem(bdys, pot)

τ    = 1e-1
sim  = NBodySimulation(sys, (0.0, 1e4τ), pbc)
sr   = run_simulation(sim, VelocityVerlet(), dt=τ)
#NBodySimulator.save_to_pdb(sr, "test.pdb")

write_xyz_traj("test.xyz", sr, 10)
