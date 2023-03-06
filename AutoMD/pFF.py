from itertools import combinations
from numpy     import exp, sqrt, dot, zeros_like, array
from numba     import njit, prange

#Morse only calculation (single molecule)
@njit(parallel=True, cache=True)
def onlyMorse(positions):
    dotSqrt = lambda diff: (sqrt(dot(diff, diff)), diff)

    #Initialize arrays and floats
    N       = len(positions)
    MorseF  = zeros_like(positions)
    ExchF   = zeros_like(positions)
    DispF   = zeros_like(positions)
    CoulF   = zeros_like(positions)
    MorseE  = ExchE = DispE = CoulE = 0.0

    #Hardcoded values
    epsilon = 11.230139012256362
    rho0    = 2.626624 # 2.3281*1.1282
    r0      = 1.1282058129221093
    preF    = 2 * epsilon * rho0 / r0

    #Iterate over all molecules
    for i in prange(0, N//2, 1):
        posi0   = positions[2*i]
        posi1   = positions[2*i+1]
        r, diff = dotSqrt(posi1 - posi0)

        #Calculate Morse Energy and Force
        expf = exp(rho0 * (1.0 - r / r0))
        E    = epsilon * expf * (expf - 2)
        F    = preF * expf * (expf - 1) * diff / r

        #Store Values
        MorseE        += E
        MorseF[2*i]   -= F
        MorseF[2*i+1] += F

    # In order to normalize the minimal morse potential to zero,
    # the following energy will be added to change the zero point
    # of all intramolecular interactions.
    MorseE += epsilon * N / 2.0

    Morse = (MorseE, MorseF)
    Exch  = ( ExchE,  ExchF)
    Disp  = ( DispE,  DispF)
    Coul  = ( CoulE,  CoulF)
    return Morse, Exch, Disp, Coul

@njit(parallel=True, cache=True)
def executeCalculations(positions, comb):
    dotSqrt = lambda diff: (sqrt(dot(diff, diff)), diff)

    def V_disp(rR, Cij):
        r2 = rR[0]*rR[0]
        E  = Cij / r2**3
        F  = 6.0 * E * rR[1] / r2
        return E, F

    def V_exch(rR, Aij, Bij):
        E = Aij * exp(-Bij * rR[0])
        F = Bij * E * rR[1] / rR[0]
        return E, F

    def calculate_Coul(rC1, rO1, rC2, rO2, rj0i0, rj1i1, rj0i1, rj1i0, ri1i0):
        V_coul  = lambda rR,Qij: ( Qij/rR[0], Qij*rR[1]/rR[0]**3)
        dotSqrt =  lambda diff: (sqrt(dot(diff, diff)), diff)

        r0          =  1.1282058129221093
        Qc          = -1.7835026375774934
        Qo          = -2.333732174702465
        alphaC      =  3.843702939952312
        alphaO      =  2.131611069944055
        wO          =  0.57135
        wC          =  0.42865

        ### CO1 ###
    	# atom positions
        r1, rCO1 = ri1i0

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
        r2, rCO2 = dotSqrt(rO2 - rC2)

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

        # rX i/j preps
        rX2i0 = dotSqrt(rX2 - rC1)
        rX2i1 = dotSqrt(rX2 - rO1)
        rj0X1 = dotSqrt(rC2 - rX1)
        rj1X1 = dotSqrt(rO2 - rX1)
        rX2X1 = dotSqrt(rX2 - rX1)

    	# C-C
        E, F    = V_coul(rj0i0, QC1*QC2)
        energy  = E
        F_Q1    = alphaC * E * rCO1 / r1
        F_Q2    = alphaC * E * rCO2 / r2
        Fi0     = - (F + F_Q1)
        Fj0     =    F - F_Q2
        Fi1     = F_Q1
        Fj1     = F_Q2

    	# C-X
        E, F    = V_coul(rX2i0, QC1*QX2)
        energy += E
        F_Q1    = alphaC * E * rCO1 / r1
        F_Q2    = - (alphaC * QC2 + alphaO * QO2) * E/QX2 * rCO2 / r2
        Fi0    -= (F + F_Q1)
        Fj0    += ((wC * F) - F_Q2)
        Fi1    += F_Q1
        Fj1    += (wO * F) + F_Q2

        # C-O
        E, F    = V_coul(rj1i0, QC1*QO2)
        energy += E
        F_Q1    = alphaC * E * rCO1 / r1
        F_Q2    = alphaO * E * rCO2 / r2
        Fi0    -= (F + F_Q1)
        Fj1    += (F + F_Q2)
        Fi1    += F_Q1
        Fj0    -= F_Q2

        # X-C
        E, F    = V_coul(rj0X1, QX1*QC2)
        energy += E
        F_Q1    = - (alphaC * QC1 + alphaO * QO1) * E/QX1 * rCO1 / r1
        F_Q2    = alphaC * E * rCO2 / r2
        Fi0    -= ((wC * F) + F_Q1)
        Fi1    += F_Q1 - (wO * F)
        Fj0    += (F - F_Q2)
        Fj1    += F_Q2

        # X-X
        E, F    = V_coul(rX2X1, QX1*QX2)
        energy += E
        F_Q1    = - (alphaC * QC1 + alphaO * QO1) * E/QX1 * rCO1 / r1
        F_Q2    = - (alphaC * QC2 + alphaO * QO2) * E/QX2 * rCO2 / r2
        Fi0    += ((-wC * F) - F_Q1)
        Fi1    += ((-wO * F) + F_Q1)
        Fj0    += (( wC * F) - F_Q2)
        Fj1    += (( wO * F) + F_Q2)


        # X-O
        E, F    = V_coul(rj1X1, QX1*QO2)
        energy += E
        F_Q1    = - (alphaC * QC1 + alphaO * QO1) * E/QX1 * rCO1 / r1
        F_Q2    = alphaO * E * rCO2 / r2
        Fi0    -= ((wC * F) + F_Q1)
        Fi1    += (F_Q1 - (wO * F))
        Fj0    -= F_Q2
        Fj1    += (F + F_Q2)

        # O-C
        E, F    = V_coul(rj0i1, QO1*QC2)
        energy += E
        F_Q1    = alphaO * E * rCO1 / r1
        F_Q2    = alphaC * E * rCO2 / r2
        Fi0    -= F_Q1
        Fi1    += (F_Q1 - F)
        Fj0    += (F - F_Q2)
        Fj1    += F_Q2

        # O-X
        E, F    = V_coul(rX2i1, QO1*QX2)
        energy += E
        F_Q1    = alphaO * E * rCO1 / r1
        F_Q2    = - (alphaC * QC2 + alphaO * QO2) * E/QX2 * rCO2 / r2
        Fi0    -= F_Q1
        Fi1    += (F_Q1 - F)
        Fj0    += ((wC * F) - F_Q2)
        Fj1    += ((wO * F) + F_Q2)

        # O-O
        E, F    = V_coul(rj1i1, QO1*QO2)
        energy += E
        F_Q1    = alphaO * E * rCO1 / r1
        F_Q2    = alphaO * E * rCO2 / r2
        Fi0    -= F_Q1
        Fi1    += (F_Q1 - F)
        Fj0    -= F_Q2
        Fj1    += (F + F_Q2)

        return energy, Fi0, Fi1, Fj0, Fj1

    #Initialize arrays and floats
    numPairs = len(comb)
    N        = len(positions)
    MorseF   = zeros_like(positions)
    ExchF    = zeros_like(positions)
    DispF    = zeros_like(positions)
    CoulF    = zeros_like(positions)
    MorseE   = ExchE = DispE = CoulE = 0.0

    #Hardcoded  Morse values
    epsilon = 11.230139012256362
    rho0    = 2.626624 # 2.3281*1.1282
    r0      = 1.1282058129221093
    preF    = 2 * epsilon * rho0 / r0

    #Hardcoded Dispersion values
    Ccc = -33.44955570065988
    Coo = -10.546349734130885
    Cco = -15.189133724736946
    Coc = -15.189133724736946

    #Hardcoded Exchange values
    Acc = 361.367206403597
    Aoo = 6370.185468304371
    Aco = 1516.76265699823
    Aoc = 1516.76265699823
    Bcc = 2.8345891887553925
    Boo = 4.2518837831330885
    Bco = 3.5432364859442407
    Boc = 3.5432364859442407

    #Iterate over all molecules
    for i in range(0, N//2, 1):
        posi0   = positions[2*i,:]
        posi1   = positions[2*i+1,:]
        r, diff = dotSqrt(posi1 - posi0)

        #Calculate Morse Energy and Force
        expf = exp(rho0 * (1.0 - r / r0))
        E    = epsilon * expf * (expf - 2)
        F    = preF * expf * (expf - 1) * diff / r

        #Store Values
        MorseE        += E
        MorseF[2*i]   -= F
        MorseF[2*i+1] += F

    # In order to normalize the minimal morse potential to zero,
    # the following energy will be added to change the zero point
    # of all intramolecular interactions.
    MorseE += epsilon * N / 2.0

    #Iterate over all pairs
    for k in prange(numPairs):
        i = comb[k][0]
        j = comb[k][1]

        posi0 = positions[2*i,:]
        posi1 = positions[2*i+1,:]
        ri1i0 = dotSqrt(posi1 - posi0)

        posj0 = positions[2*j,:]
        posj1 = positions[2*j+1,:]

        rj0i0 = dotSqrt(posj0 - posi0)
        rj1i1 = dotSqrt(posj1 - posi1)
        rj0i1 = dotSqrt(posj0 - posi1)
        rj1i0 = dotSqrt(posj1 - posi0)

        #Calculate Dispersion
        Ecc, Fcc = V_disp(rj0i0, Ccc)
        Eoo, Foo = V_disp(rj1i1, Coo)
        Eco, Fco = V_disp(rj0i1, Cco)
        Eoc, Foc = V_disp(rj1i0, Coc)

        #Store Values
        DispE                += (Ecc + Eoo + Eco + Eoc)
        DispF[2*i]           -= (Fcc + Foc)
        DispF[2*i+1]         -= (Foo + Fco)
        DispF[2*j]           += (Fcc + Fco)
        DispF[2*j+1]         += (Foo + Foc)

        #Calculate Exchange
        Ecc, Fcc  = V_exch(rj0i0, Acc, Bcc)
        Eoo, Foo  = V_exch(rj1i1, Aoo, Boo)
        Eco, Fco  = V_exch(rj0i1, Aco, Bco)
        Eoc, Foc  = V_exch(rj1i0, Aoc, Boc)

        #Store Values
        ExchE                += (Ecc + Eoo + Eco + Eoc)
        ExchF[2*i]           -= (Fcc + Foc)
        ExchF[2*i+1]         -= (Foo + Fco)
        ExchF[2*j]           += (Fcc + Fco)
        ExchF[2*j+1]         += (Foo + Foc)

        #Calculate Coulomb
        E, Fi0, Fi1, Fj0, Fj1 = calculate_Coul(posi0, posi1, posj0, posj1,
                                               rj0i0, rj1i1, rj0i1, rj1i0,
                                               ri1i0)
        #Store Values
        CoulE                += E
        CoulF[2*i]           += Fi0
        CoulF[2*i+1]         += Fi1
        CoulF[2*j]           += Fj0
        CoulF[2*j+1]         += Fj1

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
    nolabel                = True

    def __init__(self, **kwargs):
        Calculator.__init__(self, **kwargs)

    def calculate(self, atoms=None, properties=['energy'],
		  system_changes=['positions', 'numbers', 'cell',
				  'pbc', 'charges', 'magmoms']):

        Calculator.calculate(self, atoms, properties, system_changes)

        #Prep variables for doing calculations
        pos  = atoms.get_positions()

        #Get pair combinations
        comb = array([i for i in combinations(range(len(pos)//2),2)])

        if comb.size == 0:
            Morse, Exch, Disp, Coul = onlyMorse(pos)
        else:
            Morse, Exch, Disp, Coul = executeCalculations(pos, comb)

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
