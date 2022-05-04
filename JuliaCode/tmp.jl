include("./AutoMD.jl")
include("MvH_Potential.jl")
include("CustomBodies.jl")
include("CustomSimulation.jl")

using .AutoMD
using .MvH_Potential
using .CustomBodies
using .CustomSimulation
using NBodySimulator
using StaticArrays
using Unitful
using UnitfulAtomic


bdys = read_ase_xyz("10co.xyz")
pbc  = InfiniteBox()
eVA6 = 1.0u"eV*Ang^6"

M = MorseParameters(11.23u"eV", 2.3281*1.1282, 1.1282u"Ang")
d = DispersionParameters(-33.37eVA6, -10.52eVA6, -15.16eVA6)
c = CoulombParameters(3.845, 2.132, 1.1282u"Ang", false, false)
e = ExchangeParameters(361.36u"eV", 6370.10u"eV", 1516.74u"eV", 
		       2.836u"Ang^-1", 4.253u"Ang^-1", 3.544u"Ang^-1")

pot = Dict(:Morse_params      => M,
	   :Dispersion_params => d,
	   :Coulomb_params    => c,
	   :Exchange_params   => e)

sys = PotentialNBodySystem(bdys, pot)

τ   = 1e-1u"fs"
#sim  = CustomSim(sys, (0.0τ, 1e4τ))
sim = NBodySimulation(sys, (0.0τ, 1e4τ)) 
sr  = run_simulation(sim, VelocityVerlet(), dt=τ)
#NBodySimulator.save_to_pdb(sr, "test.pdb")

write_xyz_traj("test.xyz", sr, 10)
