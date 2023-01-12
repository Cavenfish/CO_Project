from os    import path
from yaml  import load
from numba import njit
from numpy import exp, sqrt, dot, zeros_like

@njit
def calculate_Morse(r, diff, D, r0, beta):
    preF = 2 * D * beta / r0
    expf = exp(beta * (r - r0))
    E    = D * (1 - 2*expf + expf**2)
    F    = preF * expf * (expf - 1) * diff / r

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

    def calculate(self, atoms=None, properties=['energy'], 
                  system_changes=all_changes):

        Calculator.calculate(self, atoms, properties, system_changes)

        #Get current atomic positions and symbols
        positions = atoms.get_positions()
        symbols   = atoms.get_chemical_symbols()

        #Execute force and energy calculations

