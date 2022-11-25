import os
os.environ["OMP_NUM_THREADS"] = '1'
os.environ["NUMBA_NUM_THREADS"] = '1'


from numpy import exp, sqrt, dot, zeros_like
from numba import njit


@njit
def calculate_Morse(r, diff):
    epsilon = 11.23
    rho0    = 2.3281*1.1282
    r0      = 1.1282
    preF    = 2 * epsilon * rho0 / r0
    expf    = exp(rho0 * (1.0 - r / r0))
    energy  = epsilon * expf * (expf - 2)
    F       = preF * expf * (expf - 1) * diff / r

    return energy, F

@njit
def V_disp(r, diff, Cij):
    r2 = r*r
    E  = Cij / r2**3
    F  = 6.0 * E * diff / r2
    return E, F

@njit
def V_exch(r, diff, Aij, Bij):
    E = Aij * exp(-Bij * r)
    F = Bij * E * diff / r
    return E, F

@njit
def calcDispEnExch(rj0i0, rj1i1, rj0i1, rj1i0):
    Ccc = -33.37
    Coo = -10.52
    Cco = -15.16
    Coc = -15.16
    Acc = 361.36
    Aoo = 6370.10
    Aco = 1516.74
    Aoc = 1516.74
    Bcc = 2.836
    Boo = 4.253
    Bco = 3.544
    Boc = 3.544

    # C-C
    DE, DFcc  = V_disp(*rj0i0, Ccc)
    EE, EFcc  = V_exch(*rj0i0, Acc, Bcc)

    # O-O
    Ed, DFoo  = V_disp(*rj1i1, Coo)
    Ee, EFoo  = V_exch(*rj1i1, Aoo, Boo)
    DE       += Ed
    EE       += Ee

    # C-O
    Ed, DFco  = V_disp(*rj0i1, Cco)
    Ee, EFco  = V_exch(*rj0i1, Aco, Bco)
    DE       += Ed
    EE       += Ee

    # O-C
    Ed, DFoc  = V_disp(*rj1i0, Coc)
    Ee, EFoc  = V_exch(*rj1i0, Aoc, Boc)
    DE       += Ed
    EE       += Ee

    #Make return tuples
    Disp = (DE, DFcc, DFoo, DFco, DFoc)
    Exch = (EE, EFcc, EFoo, EFco, EFoc)

    return Disp, Exch

@njit
def V_coul(diff, Qij):
    r = sqrt(dot(diff, diff))
    E = Qij / r
    F = Qij * diff / r**3 	# = Qij / r**2 * diff / r
    return E, F

@njit
def calculate_Coul(rC1, rO1, rC2, rO2):

    r0          =  1.1282
    Qc          = -1.786
    Qo          = -2.332
    alphaC      =  3.845
    alphaO      =  2.132
    wO          =  0.57135
    wC          =  0.42865

    ### CO1 ###
	# atom positions
    rCO1 = rO1 - rC1
    r1   = sqrt(dot(rCO1, rCO1))

    # gradient of r1 with respect to rC1: - rCO1 / r1
	# gradient of r1 with respect to rO1: rCO1 / r1

	# center of mass
    rX1 = (wC*rC1) + (wO*rO1)

	# rCO1-dependent charges:
    QC1 = Qc * exp(-alphaC*(r1 - r0))
    QO1 = Qo * exp(-alphaO*(r1 - r0))
    QX1 = - (QC1 + QO1)

    ### CO2 ###
    # atom positions
    rCO2 = rO2 - rC2
    r2   = sqrt(dot(rCO2, rCO2))

    # gradient of r2 with respect to rC2: - rCO2 / r2
    # gradient of r2 with respect to rO2: rCO2 / r2

    # center of mass
    rX2 = (wC*rC2) + (wO*rO2)

    # rCO2-dependent charges:
    QC2 = Qc * exp(-alphaC*(r2 - r0))
    QO2 = Qo * exp(-alphaO*(r2 - r0))
    QX2 = - (QC2 + QO2)

    # for force contributions resulting from bondlength-dependent charges:
    # nabla_r Q(r) = -alpha * Q(r) * nabla_r r

	# C-C
    E, F    = V_coul(rC2-rC1, QC1*QC2)
    energy  = E
    F_Q1    = alphaC * E * rCO1 / r1
    F_Q2    = alphaC * E * rCO2 / r2
    Fi0     = - (F + F_Q1)
    Fj0     =    F - F_Q2
    Fi1     = F_Q1
    Fj1     = F_Q2

	# C-X
    E, F    = V_coul(rX2-rC1, QC1*QX2)
    energy += E
    F_Q1    = alphaC * E * rCO1 / r1
    F_Q2    = - (alphaC * QC2 + alphaO * QO2) * E/QX2 * rCO2 / r2
    Fi0    -= (F + F_Q1)
    Fj0    += ((wC * F) - F_Q2)
    Fi1    += F_Q1
    Fj1    += (wO * F) + F_Q2

    # C-O
    E, F    = V_coul(rO2-rC1, QC1*QO2)
    energy += E
    F_Q1    = alphaC * E * rCO1 / r1
    F_Q2    = alphaO * E * rCO2 / r2
    Fi0    -= (F + F_Q1)
    Fj1    += (F + F_Q2)
    Fi1    += F_Q1
    Fj0    -= F_Q2

    # X-C
    E, F    = V_coul(rC2-rX1, QX1*QC2)
    energy += E
    F_Q1    = - (alphaC * QC1 + alphaO * QO1) * E/QX1 * rCO1 / r1
    F_Q2    = alphaC * E * rCO2 / r2
    Fi0    -= ((wC * F) + F_Q1)
    Fi1    += F_Q1 - (wO * F)
    Fj0    += (F - F_Q2)
    Fj1    += F_Q2

    # X-X
    E, F    = V_coul(rX2-rX1, QX1*QX2)
    energy += E
    F_Q1    = - (alphaC * QC1 + alphaO * QO1) * E/QX1 * rCO1 / r1
    F_Q2    = - (alphaC * QC2 + alphaO * QO2) * E/QX2 * rCO2 / r2
    Fi0    += ((-wC * F) - F_Q1)
    Fi1    += ((-wO * F) + F_Q1)
    Fj0    += (( wC * F) - F_Q2)
    Fj1    += (( wO * F) + F_Q2)


    # X-O
    E, F    = V_coul(rO2-rX1, QX1*QO2)
    energy += E
    F_Q1    = - (alphaC * QC1 + alphaO * QO1) * E/QX1 * rCO1 / r1
    F_Q2    = alphaO * E * rCO2 / r2
    Fi0    -= ((wC * F) + F_Q1)
    Fi1    += (F_Q1 - (wO * F))
    Fj0    -= F_Q2
    Fj1    += (F + F_Q2)

    # O-C
    E, F    = V_coul(rC2-rO1, QO1*QC2)
    energy += E
    F_Q1    = alphaO * E * rCO1 / r1
    F_Q2    = alphaC * E * rCO2 / r2
    Fi0    -= F_Q1
    Fi1    += (F_Q1 - F)
    Fj0    += (F - F_Q2)
    Fj1    += F_Q2

    # O-X
    E, F    = V_coul(rX2-rO1, QO1*QX2)
    energy += E
    F_Q1    = alphaO * E * rCO1 / r1
    F_Q2    = - (alphaC * QC2 + alphaO * QO2) * E/QX2 * rCO2 / r2
    Fi0    -= F_Q1
    Fi1    += (F_Q1 - F)
    Fj0    += ((wC * F) - F_Q2)
    Fj1    += ((wO * F) + F_Q2)

    # O-O
    E, F    = V_coul(rO2-rO1, QO1*QO2)
    energy += E
    F_Q1    = alphaO * E * rCO1 / r1
    F_Q2    = alphaO * E * rCO2 / r2
    Fi0    -= F_Q1
    Fi1    += (F_Q1 - F)
    Fj0    -= F_Q2
    Fj1    += (F + F_Q2)

    return energy, Fi0, Fi1, Fj0, Fj1

@njit
def diffDotSqrt(v2, v1):
    diff = v2 - v1
    r    = sqrt(dot(diff, diff))
    return (r, diff)

@njit(parallel=False)
def executeCalculations(positions):
    N      = len(positions)
    MorseF = ExchF = DispF = CoulF = zeros_like(positions)
    MorseE = ExchE = DispE = CoulE = 0.0

    for i in range(0, N//2, 1):
        posi0 = positions[2*i]
        posi1 = positions[2*i+1]
        ri1i0 = diffDotSqrt(posi1, posi0)

        #Calculate Morse
        E, F           = calculate_Morse(*ri1i0)
        MorseE        += E
        MorseF[2*i]   -= F
        MorseF[2*i+1] += F

        for j in range(i+1, N//2, 1):
            posj0 = positions[2*j]
            posj1 = positions[2*j+1]

            rj0i0 = diffDotSqrt(posj0, posi0)
            rj1i1 = diffDotSqrt(posj1, posi1)
            rj0i1 = diffDotSqrt(posj0, posi1)
            rj1i0 = diffDotSqrt(posj1, posi0)

            #Calculate Dispersion and Exchange
            disp, exch            = calcDispEnExch(rj0i0, rj1i1, rj0i1, rj1i0)

            #Store Dispersion
            E, Fcc, Foo, Fco, Foc = disp
            DispE                += E
            DispF[2*i]           -= (Fcc + Foc)
            DispF[2*i+1]         -= (Foo + Fco)
            DispF[2*j]           += (Fcc + Fco)
            DispF[2*j+1]         += (Foo + Foc)

            #Store Exchange
            E, Fcc, Foo, Fco, Foc = exch
            ExchE                += E
            ExchF[2*i]           -= (Fcc + Foc)
            ExchF[2*i+1]         -= (Foo + Fco)
            ExchF[2*j]           += (Fcc + Fco)
            ExchF[2*j+1]         += (Foo + Foc)

            #Calculate Coulomb
            E, Fi0, Fi1, Fj0, Fj1 = calculate_Coul(posi0, posi1, posj0, posj1)
            CoulE                += E
            CoulF[2*i]           += Fi0
            CoulF[2*i+1]         += Fi1
            CoulF[2*j]           += Fj0
            CoulF[2*j+1]         += Fj1

    # In order to normalize the minimal morse potential to zero,
    # the following energy will be added to change the zero point
    # of all intramolecular interactions.
    MorseE += 11.230 * N / 2.0

    Morse = (MorseE, MorseF)
    Exch  = ( ExchE,  ExchF)
    Disp  = ( DispE,  DispF)
    Coul  = ( CoulE,  CoulF)
    return Morse, Exch, Disp, Coul


#from ase.units import Bohr,Hartree
from ase.calculators.calculator import Calculator

class MvH_CO(Calculator):
    """
    CO-CO interaction potential as described in:
	Karssemeijer, L. J.; Ioppolo, S.; van Hemert, M. C.;
        van der Avoird, A.; Allodi, M. A.; Blake, G. A.; Cuppen, H. M.
        Dynamics of CO in Amorphous Water-Ice Environments.
        ApJ 2013, 781 (1), 16. https://doi.org/10.1088/0004-637X/781/1/16.

    or likewise (in atomic units)
	van Hemert, M. C.; Takahashi, J.; van Dishoeck, E. F.
	Molecular Dynamics Study of the Photodesorption of CO Ice.
	J. Phys. Chem. A 2015, 119 (24), 6354–6369.
	https://doi.org/10.1021/acs.jpca.5b02611.
    """

    implemented_properties = ['energy', 'forces']

    def __init__(self, **kwargs):
        Calculator.__init__(self, **kwargs)

    def calculate(self, atoms=None, properties=['energy'],
		  system_changes=['positions', 'numbers', 'cell',
				  'pbc', 'charges', 'magmoms']):

        Calculator.calculate(self, atoms, properties, system_changes)

        #Prep variables for doing calculations
        positions    = atoms.get_positions()

        #Execute calculations
        Morse, Exch, Disp, Coul = executeCalculations(positions)

        E_pair = Exch[0] + Disp[0] + Coul[0]
        F_pair = Exch[1] + Disp[1] + Coul[1]

        self.results['energy']   = Morse[0] + E_pair
        self.results['forces']   = Morse[1] + F_pair

        self.results['E_Morse'] = Morse[0]
        self.results['E_Disp']  =  Disp[0]
        self.results['E_Exch']  =  Exch[0]
        self.results['E_Coul']	=  Coul[0]
        self.results['E_pair']  = E_pair

        self.results['F_Morse'] = Morse[1]
        self.results['F_Disp']  =  Disp[1]
        self.results['F_Exch']  =  Exch[1]
        self.results['F_Coul']	=  Coul[1]
        self.results['F_pair']  = F_pair

    def get_energy_contributions(self):
        E_intra  = self.results['E_Morse']
        E_pair   = self.results['E_pair']
        E_ex     = self.results['E_Exch']
        E_disp   = self.results['E_Disp']
        E_elst   = self.results['E_Coul']
        return E_intra, E_pair, E_ex, E_disp, E_elst

    def get_force_contributions(self):
        F_intra  = self.results['F_Morse']
        F_pair   = self.results['F_pair']
        F_ex     = self.results['F_Exch']
        F_disp   = self.results['F_Disp']
        F_elst   = self.results['F_Coul']
        return F_intra, F_pair, F_ex, F_disp, F_elst
