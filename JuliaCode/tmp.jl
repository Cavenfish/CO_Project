include("./AutoMD.jl")
include("MvH_Potential.jl")
include("CustomBodies.jl")

using .AutoMD
using .MvH_Potential
using .CustomBodies
using NBodySimulator
using StaticArrays
using Unitful
using UnitfulAtomic

# Time to beat for 100fs run => 33seconds 

@time begin
bdys = read_ase_xyz("250co0.xyz")
pbc  = InfiniteBox()
end

M = MorseParameters(11.23, 2.3281*1.1282, 1.1282)
d = DispersionParameters(-33.37, -10.52, -15.16)
c = CoulombParameters(3.845, 2.132, 1.1282, true, true)
e = ExchangeParameters(361.36, 6370.10, 1516.74, 2.836, 4.253, 3.544)

pot = Dict(:Morse_params      => M,
	   :Dispersion_params => d,
	   :Coulomb_params    => c,
	   :Exchange_params   => e)

sys = PotentialNBodySystem(bdys, pot)

τ   = uconvert(u"Ang*u^0.5*eV^-0.5", 1u"fs")
τ   = ustrip(τ)
sim = NBodySimulation(sys, (0.0τ, 1e2τ)) 
@time begin
sr  = run_simulation(sim, VelocityVerlet(), dt=τ)
end
#NBodySimulator.save_to_pdb(sr, "test.pdb")
@time begin
write_xyz_traj("test.xyz", sr, 10)
end
