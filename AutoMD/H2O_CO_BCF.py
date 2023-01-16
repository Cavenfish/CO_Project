from os    import path
from yaml  import load
from numba import njit
from numpy import exp, sqrt, dot, zeros_like, where

@njit
def calculate_Morse(r, diff, D, r0, beta):
    preF = 2 * D * beta / r0
    expf = exp(beta * (r - r0))
    E    = D * (1 - 2*expf + expf**2)
    F    = preF * expf * (expf - 1) * diff / r
    #TODO: normalize values here 
    return E, F

@njit
def V_HOH(K, theta, theta0):
    E = 0.5 * K * (theta - theta0)**2
    F = 0
    return E, F

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
def V_coul(r, diff, Qij):
    E = Qij / r
    F = Qij * diff / r**3
    return E, F

@njit
def V_LJ(r, diff, sigma, epsilon):
    E = 4 * epsilon * ( (sigma/r)**12 - (sigma/r)**6 )
    c = (4 * epsilon * sigma) / r**3
    k = 12 * (sigma/r)**11 - 6 * (sigma/r)**5
    F = c * k * diff
    return E, F

@njit
def calculate_COCO(molA, molB):
    return

@njit
def calculate_H2OH2O(molA, molB):
    return

@njit
def calculate_COH2O(molA, molB):
    return

@njit
def calculate_CO(pos, molA, params, MorseE, MorseF):
    posi0 = positions[molA[0]]
    posi1 = positions[molA[1]]
    ri1i0 = diffDotSqrt(posi1, posi0)

    intraP = params['CO'] if molA in co else params['H2O']
    Dr0B   = (intraP['D'], intraP['r0'], intraP['beta'])

    #Calculate Morse 
    E, F             = calculate_Morse(*ri1i0, *Dr0B)
    MorseE          += E
    MorseF[molA[0]] -= F
    MorseF[molA[1]] += F

    return MorseE, MorseF

@njit
def diffDotSqrt(v2, v1):
    diff = v2 - v1
    r    = sqrt(dot(diff,diff))
    return (r, diff)

@njit
def executeCalculations(positions, co, water, params):
    molecs = co.append(water)
    N      = len(molecs)
    MorseF = zeros_like(positions)
    HarmoF = zeros_like(positions)
    ExchF  = zeros_like(positions)
    DispF  = zeros_like(positions)
    CoulF  = zeros_like(positions)
    MorseE = HarmoE = ExchE = DispE = CoulE = 0.0
    
    for i in range(N):
        molA  = molecules[i]
        
        if molA in co:
            MorseE, MorseF = calculate_CO(positions, molA, params,
                                          MorseE, MorseF)
        else:
            MorseE, = calculate_H2O(positions, molA, params)



    return

from ase.calculators.calculator import (Calculator, all_changes)

class H2O_CO(Calculator):
    """
    H2O-CO interaction potential as described in:

    """

    implemented_properties = ['energy', 'forces']
    nolabel                = True
    
    #Load FF params
    fname = path.dirname(__file__) + '/ff_params.yaml'
    with open(fname, 'r') as f:
        params = load(f)

    def __init__(self, **kwargs):
        Calculator.__init__(self, **kwargs)

    def get_molecules(self, atoms):
        def get_co(i,d):
            j = where((d[i] > 0.1) & (d[i] < 1.3))[0][0]
            return (i,j)

        def get_water(i,d):
            x = where((d[i] > 0.1) & (d[i] < 1.9))[0] 
            m = [i,x[0],x[1]]
            m.sort()
            return tuple(m)

        symbols   = atoms.get_chemical_symbols()
        distances = atoms.get_all_distances()
        N         = range(len(symbols))
        carbons   = [i for i in N if symbols[i] == 'C']
        hydrogen  = [i for i in N if symbols[i] == 'H']
        co        = [get_co(i,distances) for i in carbons]
        water     = {get_water(i,distances) for i in hydrogen}

        return co, water

    def calculate(self, atoms=None, properties=['energy'], 
                  system_changes=all_changes):

        Calculator.calculate(self, atoms, properties, system_changes)

        #Get current atomic positions and symbols
        positions = atoms.get_positions()
        symbols   = atoms.get_chemical_symbols()

        #Execute force and energy calculations

