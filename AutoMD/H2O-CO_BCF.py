import yaml
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
    F =
    return E, F

@njit
def V_coul(r, diff, Qij):
    E = Qij / r
    F = Qij * diff / r**3
    return E, F

@njit
def V_LJ(r, diff, sigma, epsilon):
    E = 4 * epsilon * ( (sigma/r)**12 - (sigma/r)**6 )
    F = 
    return E, F
