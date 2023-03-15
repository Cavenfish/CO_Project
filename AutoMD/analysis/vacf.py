from ..config import *
from scipy import fft
from numpy.fft import fftfreq
from numba import njit

@njit
def vacf(vel, m):
    natoms = vel.shape[1]
    nstep  = vel.shape[0]
    c      = np.zeros(nstep)
    for t in range(nstep):
        for j in range(nstep-t):
            for i in range(natoms):
                c[t] += np.dot(vel[j,i,:], vel[j+t,i,:]) * m[i]
    
    c /= c[0]
    return c

def vdos(c, dt):
    k = 29979245800.0 # Divide by this to convert from Hz to cm^-1
    I = fft(c)
    v = fftfreq(len(c), dt) / k
    I = I[:len(v)//2]
    v = v[:len(v)//2]
    return v,I

def window():
    return
