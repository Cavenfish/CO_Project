import numpy as np
from numpy import exp, sqrt

def calculate_Morse(positions, epsilon, rho0, r0):
    energy = 0.0
    forces = np.zeros_like(positions)
    preF   = 2 * epsilon * rho0 / r0
    N      = len(positions)
    E_molc = np.zeros(N//2)

    for i in range(0, N//2, 1):
        diff           = positions[2*i+1] - positions[2*i]
        r              = sqrt(np.dot(diff, diff))
        expf           = exp(rho0 * (1.0 - r / r0))
        energy        += epsilon * expf * (expf - 2)
        F              = preF * expf * (expf - 1) * diff / r
        forces[2*i]   -= F
        forces[2*i+1] += F
        E_molc[i]     += (epsilon * expf * (expf - 2)) + epsilon
    
    # In order to normalize the minimal morse potential to zero,
    # the following energy will be added to change the zero point
    # of all intramolecular interactions.
    energy += epsilon * N / 2.0
    return energy, forces, E_molc

def V_disp(diff, Cij):
    r2 = np.dot(diff, diff)
    E  = Cij / r2**3
    F  = 6.0 * E * diff / r2
    return E, F

def calculate_Disp(positions, Ccc,Coo,Cco,Coc):
    Ecc = Eoo = Eco = Eoc = 0.0
    forces = np.zeros_like(positions)
    N      = len(positions)
    E_molc = np.zeros(N//2)

    for i in range(0,N//2,1):
        for j in range(i+1,N//2,1):
            # C-C
            E, F         = V_disp(positions[2*j] - positions[2*i], Ccc)
            Ecc         += E
            forces[2*i] -= F
            forces[2*j] += F
            E_molc[i]   += E
	    
            # O-O
            E, F           = V_disp(positions[2*j+1] - positions[2*i+1], Coo)
            Eoo           += E
            forces[2*i+1] -= F
            forces[2*j+1] += F
            E_molc[i]     += E
	    
            # C-O
            E, F           = V_disp(positions[2*j] - positions[2*i+1], Cco)
            Eco           += E
            forces[2*i+1] -= F
            forces[2*j]   += F
            E_molc[i]     += E
	    
            # O-C
            E, F           = V_disp(positions[2*j+1] - positions[2*i], Coc)
            Eoc           += E
            forces[2*i]   -= F
            forces[2*j+1] += F
            E_molc[i]     += E
	
    energy = Ecc + Eoo + Eco + Eoc
    return energy, forces, E_molc

def V_exch(diff, Aij, Bij):
    r = sqrt(np.dot(diff, diff))
    E = Aij * exp(-Bij * r)
    F = Bij * E * diff / r
    return E, F

def calculate_Exch(positions, Acc,Aoo,Aco,Aoc, Bcc,Boo,Bco,Boc):
    Ecc = Eoo = Eco = Eoc = 0.0
    forces = np.zeros_like(positions)
    N      = len(positions)
    E_molc = np.zeros(N//2)
    
    for i in range(0,N//2,1):
        for j in range(i+1,N//2,1):
            
            # C-C
            E, F         = V_exch(positions[2*j] - positions[2*i], Acc, Bcc)
            Ecc         += E
            forces[2*i] -= F
            forces[2*j] += F
            E_molc[i]   += E
	    
            # O-O
            E, F           = V_exch(positions[2*j+1] - positions[2*i+1], 
                                    Aoo, Boo)
            Eoo           += E
            forces[2*i+1] -= F
            forces[2*j+1] += F
            E_molc[i]     += E
			
            # C-O
            E, F           = V_exch(positions[2*j] - positions[2*i+1], Aco, Bco)
            Eco           += E
            forces[2*i+1] -= F
            forces[2*j]   += F
            E_molc[i]     += E
	    
            # O-C
            E, F           = V_exch(positions[2*j+1] - positions[2*i], Aoc, Boc)
            Eoc           += E
            forces[2*i]   -= F
            forces[2*j+1] += F
            E_molc[i]     += E
    
    energy = Ecc + Eoo + Eco + Eoc
    return energy, forces, E_molc

def V_coul(diff, Qij):
    r = sqrt(np.dot(diff, diff))
    E = Qij / r
    F = Qij * diff / r**3 	# = Qij / r**2 * diff / r
    return E, F

def calculate_Coul(positions, cms_weights, r0,Qc,Qo,alphaC,alphaO, coco_FQ_mode=True, coco_FX_mode=True):

    Ecc = Eoo = Eco = Eoc = 0.0
    Ecx = Exc = Eox = Exo = Exx = 0.0
    forces = np.zeros_like(positions)
    N      = len(positions)
    E_molc = np.zeros(N//2)

    forces_c_o   = np.zeros_like(positions)
    forces_c_o_Q = np.zeros_like(positions)
    forces_x     = np.zeros((N//2, 3))
    forces_x_Q   = np.zeros_like(forces_x)

    for i in range(0,N//2,1):

        ### CO1 ###
	# atom positions
        rC1  = positions[2*i]
        rO1  = positions[2*i+1]
        rCO1 = rO1 - rC1
        r1   = sqrt(np.dot(rCO1, rCO1))
	
        # gradient of r1 with respect to rC1: - rCO1 / r1
	# gradient of r1 with respect to rO1: rCO1 / r1

	# center of mass
        wC1 = cms_weights[2*i]
        wO1 = cms_weights[2*i+1]
        rX1 = (wC1*rC1) + (wO1*rO1)
	#wFX1 = 1.0/wO1 - 1.0/wC1

	# rCO1-dependent charges:
        QC1 = Qc * exp(-alphaC*(r1 - r0))
        QO1 = Qo * exp(-alphaO*(r1 - r0))
        QX1 = - (QC1 + QO1)

        for j in range(i+1,N//2,1):

            ### CO2 ###
	    # atom positions
            rC2  = positions[2*j]
            rO2  = positions[2*j+1]
            rCO2 = rO2 - rC2
            r2   = sqrt(np.dot(rCO2, rCO2))
	    
            # gradient of r2 with respect to rC2: - rCO2 / r2
	    # gradient of r2 with respect to rO2: rCO2 / r2

	    # center of mass
            wC2 = cms_weights[2*j]
            wO2 = cms_weights[2*j+1]
            rX2 = (wC2*rC2) + (wO2*rO2)
	    #wFX2 = 1.0/wO2 - 1.0/wC2

	    # rCO2-dependent charges:
            QC2 = Qc * exp(-alphaC*(r2 - r0))
            QO2 = Qo * exp(-alphaO*(r2 - r0))
            QX2 = - (QC2 + QO2)

	    # for force contributions resulting from bondlength-dependent charges:
	    # nabla_r Q(r) = -alpha * Q(r) * nabla_r r

	    # C-C
            E, F                 = V_coul(rC2-rC1, QC1*QC2)
            Ecc                 += E
            forces_c_o[2*i]     -= F
            forces_c_o[2*j]     += F
            F_Q1                 = alphaC * E * rCO1 / r1
            F_Q2                 = alphaC * E * rCO2 / r2
            forces_c_o_Q[2*i]   -= F_Q1
            forces_c_o_Q[2*i+1] += F_Q1
            forces_c_o_Q[2*j]   -= F_Q2
            forces_c_o_Q[2*j+1] += F_Q2
            E_molc[i]           += E

	    # C-X
            E, F                 = V_coul(rX2-rC1, QC1*QX2)
            Ecx                 += E
            forces_c_o[2*i]     -= F
            forces_x[j]         += F
            forces_c_o[2*j]     += wC2 * F
            forces_c_o[2*j+1]   += wO2 * F
            F_Q1                 = alphaC * E * rCO1 / r1
            F_Q2                 = - (alphaC * QC2 + alphaO * QO2) * E/QX2 * rCO2 / r2
            forces_c_o_Q[2*i]   -= F_Q1
            forces_c_o_Q[2*i+1] += F_Q1
	    #forces_x_Q[j]       += wFX2 * F_Q2
            forces_c_o_Q[2*j]   -= F_Q2
            forces_c_o_Q[2*j+1] += F_Q2
            E_molc[i]           += E

	    # C-O
            E, F                 = V_coul(rO2-rC1, QC1*QO2)
            Eco                 += E
            forces_c_o[2*i]     -= F
            forces_c_o[2*j+1]   += F
            F_Q1                 = alphaC * E * rCO1 / r1
            F_Q2                 = alphaO * E * rCO2 / r2
            forces_c_o_Q[2*i]   -= F_Q1
            forces_c_o_Q[2*i+1] += F_Q1
            forces_c_o_Q[2*j]   -= F_Q2
            forces_c_o_Q[2*j+1] += F_Q2
            E_molc[i]           += E

	    # X-C
            E, F                 = V_coul(rC2-rX1, QX1*QC2)
            Exc                 += E
            forces_x[i]         -= F
            forces_c_o[2*i]     -= wC1 * F
            forces_c_o[2*i+1]   -= wO1 * F
            forces_c_o[2*j]     += F
            F_Q1                 = - (alphaC * QC1 + alphaO * QO1) * E/QX1 * rCO1 / r1
            F_Q2                 = alphaC * E * rCO2 / r2
	    #forces_x_Q[i]       += wFX1 * F_Q1
            forces_c_o_Q[2*i]   -= F_Q1
            forces_c_o_Q[2*i+1] += F_Q1
            forces_c_o_Q[2*j]   -= F_Q2
            forces_c_o_Q[2*j+1] += F_Q2
            E_molc[i]           += E

	    # X-X
            E, F                 = V_coul(rX2-rX1, QX1*QX2)
            Exx                 += E
            forces_x[i]         += -F
            forces_c_o[2*i]     += -wC1 * F
            forces_c_o[2*i+1]   += -wO1 * F
            forces_x[j]         +=  F
            forces_c_o[2*j]     += wC2 * F
            forces_c_o[2*j+1]   += wO2 * F
            F_Q1                 = - (alphaC * QC1 + alphaO * QO1) * E/QX1 * rCO1 / r1
            F_Q2                 = - (alphaC * QC2 + alphaO * QO2) * E/QX2 * rCO2 / r2
	    #forces_x[i]         += wFX1 * F_Q1
            forces_c_o_Q[2*i]   -= F_Q1
            forces_c_o_Q[2*i+1] += F_Q1
	    #forces_x[j]         += wFX2 * F_Q2
            forces_c_o_Q[2*j]   -= F_Q2
            forces_c_o_Q[2*j+1] += F_Q2
            E_molc[i]           += E

	    # X-O
            E, F                 = V_coul(rO2-rX1, QX1*QO2)
            Exo                 += E
            forces_x[i]         -= F
            forces_c_o[2*i]     -= wC1 * F
            forces_c_o[2*i+1]   -= wO1 * F
            forces_c_o[2*j+1]   += F
            F_Q1                 = - (alphaC * QC1 + alphaO * QO1) * E/QX1 * rCO1 / r1
            F_Q2                 = alphaO * E * rCO2 / r2
	    #forces_x[i]        += wFX1 * F_Q1
            forces_c_o_Q[2*i]   -= F_Q1
            forces_c_o_Q[2*i+1] += F_Q1
            forces_c_o_Q[2*j]   -= F_Q2
            forces_c_o_Q[2*j+1] += F_Q2
            E_molc[i]           += E

	    # O-C
            E, F                 = V_coul(rC2-rO1, QO1*QC2)
            Eoc                 += E
            forces_c_o[2*i+1]   -= F
            forces_c_o[2*j]     += F
            F_Q1                 = alphaO * E * rCO1 / r1
            F_Q2                 = alphaC * E * rCO2 / r2
            forces_c_o_Q[2*i]   -= F_Q1
            forces_c_o_Q[2*i+1] += F_Q1
            forces_c_o_Q[2*j]   -= F_Q2
            forces_c_o_Q[2*j+1] += F_Q2
            E_molc[i]           += E

	    # O-X
            E, F                 = V_coul(rX2-rO1, QO1*QX2)
            Eox                 += E
            forces_c_o[2*i+1]   -= F
            forces_x[j]         += F
            forces_c_o[2*j]     += wC2 * F
            forces_c_o[2*j+1]   += wO2 * F
            F_Q1                 = alphaO * E * rCO1 / r1
            F_Q2                 = - (alphaC * QC2 + alphaO * QO2) * E/QX2 * rCO2 / r2
            forces_c_o_Q[2*i]   -= F_Q1
            forces_c_o_Q[2*i+1] += F_Q1
            #forces_x_Q[j]       += wFX2 * F_Q2
            forces_c_o_Q[2*j]   -= F_Q2
            forces_c_o_Q[2*j+1] += F_Q2
            E_molc[i]           += E

	    # O-O
            E, F                 = V_coul(rO2-rO1, QO1*QO2)
            Eoo                 += E
            forces_c_o[2*i+1]   -= F
            forces_c_o[2*j+1]   += F
            F_Q1                 = alphaO * E * rCO1 / r1
            F_Q2                 = alphaO * E * rCO2 / r2
            forces_c_o_Q[2*i]   -= F_Q1
            forces_c_o_Q[2*i+1] += F_Q1
            forces_c_o_Q[2*j]   -= F_Q2
            forces_c_o_Q[2*j+1] += F_Q2
            E_molc[i]           += E

    energy = Ecc + Eoo + Eco + Eoc + \
			 Ecx + Exc + Eox + Exo + Exx

    if coco_FQ_mode:
        forces = forces_c_o
    else:
        forces = forces_c_o + forces_c_o_Q

    if not coco_FX_mode:
	# Subtract forces parallel to bonds from X forces
        forces_to_redistribute = np.zeros_like(forces_x)
        for i in range(0, N//2, 1):
            rC1                       = positions[2*i]
            rO1                       = positions[2*i+1]
            rCO                       = rC1 - rO1
            eCO                       = rCO / np.linalg.norm(rCO)
            FX                        = forces_x[i] + forces_x_Q[i]
            FX_parallel               = np.dot(FX,eCO) * eCO
            FX_perp                   = FX - FX_parallel
            forces_to_redistribute[i] = FX_perp
	    
            # redistribute X forces
            for i in range(0, N//2, 1):
                wC               = cms_weights[2*i]
                wO               = cms_weights[2*i+1]
                forces[2*i,:]   += wC * forces_to_redistribute[i]
                forces[2*i+1,:] += wO * forces_to_redistribute[i]

    return energy, forces, E_molc


from ase.units import Bohr,Hartree
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

    implemented_properties = ['energy', 'forces', 'energies']

    # epsilon: Absolute minimum depth, default 11.230
    # rho0   : Exponential prefactor
    # r0     : Minimum distance, default 1.1282
    # Cij    : Dispersion coefficients in eV*Angstrom^6
    parameters_FS = {   'epsilon':   11.230,         
			'rho0'   :    2.3281*1.1282,
			'r0'     :    1.1282,
			'Ccc'    :  -33.37,
			'Coo'    :  -10.52,
			'Cco'    :  -15.16,
			'Coc'    :  -15.16,
			'Acc'    :  361.36,
			'Aoo'    : 6370.10,
			'Aco'    : 1516.74,
			'Aoc'    : 1516.74,
			'Bcc'    :    2.836,
			'Boo'    :    4.253,
			'Bco'    :    3.544,
			'Boc'    :    3.544,
		        'Qc'     :   -1.786,
			'Qo'     :   -2.332,
			'alphaC' :    3.845,
		        'alphaO' :    2.132 }
    
    # The force constant in the potential minimum is k = 2 * epsilon * (rho0 / r0)**2, default 2.3281*1.1282

    # taken from cocdyn/blkdtexp.f
    rcoeq  =   2.132
    #fcco  =  12.2	# harmonic force constant
    facc   =   0.57135
    faco   =   0.42865
    beta   =   1.232
    edis   =   0.4127
    
    # charges
    qec    =  -0.47
    qeo    =  -0.615	
    alphac =   2.034
    alphao =   1.128
    
    # CC
    c6cc   = -55.98
    acc    =  13.28
    alcc   =   1.50
    
    # OO
    c6oo   = -17.65
    aoo    = 234.1
    aloo   =   2.250
	
    # OC
    c6oc   = -25.42
    aoc    =  55.74
    aloc   =   1.875
    
    # CO
    c6co   = -25.42
    aco    =  55.74
    alco   =   1.875
    
    parameters_coco = { 'epsilon': edis   * Hartree,
			'rho0'   : beta   * rcoeq,
			'r0'     : rcoeq  * Bohr,
			'Ccc'    : c6cc   * Hartree * Bohr**6,
			'Coo'    : c6oo   * Hartree * Bohr**6,
			'Cco'    : c6co   * Hartree * Bohr**6,
			'Coc'    : c6oc   * Hartree * Bohr**6,
			'Acc'    : acc    * Hartree,
			'Aoo'    : aoo    * Hartree,
			'Aco'    : aco    * Hartree,
			'Aoc'    : aoc    * Hartree,
			'Bcc'    : alcc   / Bohr,
			'Boo'    : aloo   / Bohr,
			'Bco'    : alco   / Bohr,
			'Boc'    : aloc   / Bohr,
			'Qc'     : qec    * Hartree**0.5 * Bohr**0.5,
			'Qo'     : qeo    * Hartree**0.5 * Bohr**0.5,
			'alphaC' : alphac / Bohr,
			'alphaO' : alphao / Bohr}

    # COCO_Q_VARIABLE  ->  If True, use variable charges
    # COCO_CMS_WEIGHTS ->  If True, use center-of-mass weights from coco
    # COCO_FQ_MODE     ->  If True, do not include force contributions due
    #                        to CO-bond-length-dependent charges
    # COCO_FX_MODE     ->  If True, do not add torque-free projection of
    #                        forces on X to C and O 
    parameters_runtime = { "COCO_Q_VARIABLE"  : True,
	                   "COCO_CMS_WEIGHTS" : True,
	                   "COCO_FQ_MODE"     : False,
	                   "COCO_FX_MODE"     : True}

    default_parameters = dict(parameters_coco, **parameters_runtime)
    nolabel            = True

    def __init__(self, **kwargs):
        Calculator.__init__(self, **kwargs)

    def get_cms_weights(self):
        masses = self.atoms.get_masses()
        N      = len(masses)
	
        # ignoring ASE masses - not sensitive to isotopes
        if self.parameters.COCO_CMS_WEIGHTS:
            cms_weights = np.array( [self.faco, self.facc]*N )
        else:
            cms_weights = np.zeros_like(masses)
            for i in range(0, N, 2):
                mC, mO           = masses[i:i+2]
                MCO              = mC + mO
                cms_weights[i]   = mC / MCO
                cms_weights[i+1] = mO / MCO
	
        return cms_weights


    def calculate(self, atoms=None, properties=['energy'],
		  system_changes=['positions', 'numbers', 'cell',
				  'pbc', 'charges', 'magmoms']):

        Calculator.calculate(self, atoms, properties, system_changes)

        #Prep variables for doing calculations
        positions    = atoms.get_positions()
        morse_list   = ('epsilon', 'rho0', 'r0')
        exch_list    = ('Acc', 'Aoo', 'Aco', 'Aoc', 'Bcc', 'Boo', 'Bco', 'Boc')
        disp_list    = ('Ccc', 'Coo', 'Cco', 'Coc')
        Morse_kwargs = {k:self.parameters[k] for k in morse_list}
        Exch_kwargs  = {k:self.parameters[k] for k in  exch_list}
        Disp_kwargs  = {k:self.parameters[k] for k in  disp_list}

        #Execute calculations
        E_Morse, F_Morse, E_molc_Morse = calculate_Morse(positions,
                                                         **Morse_kwargs)
        E_Exch, F_Exch, E_molc_Exch    = calculate_Exch( positions, 
                                                         **Exch_kwargs)
        E_Disp, F_Disp, E_molc_Disp    = calculate_Disp( positions, 
                                                         **Disp_kwargs)

        cms_weights = self.get_cms_weights()
        coul_list   = ('r0','Qc','Qo','alphaC','alphaO')
        Coul_kwargs = {k:self.parameters[k] for k in coul_list}

        if not self.parameters.COCO_Q_VARIABLE:
            Coul_kwargs['alphaC']   = 0.0
            Coul_kwargs['alphaO']   = 0.0
        Coul_kwargs['coco_FQ_mode'] = self.parameters.COCO_FQ_MODE
        Coul_kwargs['coco_FX_mode'] = self.parameters.COCO_FX_MODE
        E_Coul, F_Coul, E_molc_Coul = calculate_Coul(positions, 
                                                     cms_weights, 
                                                     **Coul_kwargs)

        E_pair = E_Exch + E_Disp + E_Coul
        F_pair = F_Exch + F_Disp + F_Coul
        E_molc = E_molc_Morse + E_molc_Exch + E_molc_Disp + E_molc_Coul

        self.results['energy']   = E_Morse + E_pair
        self.results['forces']   = F_Morse + F_pair
        self.results['energies'] = E_molc

        self.results['E_Morse'] = E_Morse
        self.results['E_Disp']  = E_Disp
        self.results['E_Exch']  = E_Exch
        self.results['E_Coul']	= E_Coul
        self.results['E_pair']  = E_pair

        self.results['F_Morse'] = F_Morse
        self.results['F_Disp']  = F_Disp
        self.results['F_Exch']  = F_Exch
        self.results['F_Coul']	= F_Coul
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
