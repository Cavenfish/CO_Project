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
                c[t] += np.dot(vel[j,i], vel[j+t,i]) * m[i]

    c /= c[0]
    return c

class VACF:
    '''
    Velocity AutoCorrelation Function
    '''
    #Window Functions
    Hann       = lambda _,n,N: np.sin((np.pi*n)/N)**2
    Welch      = lambda _,n,N: 1 - ((n-N/2)/(N/2))**2
    HannMirror = lambda _,n,N: np.cos((np.pi*n)/(2*(N-1)))**2

    
    def  __init__(self, vel, dt, m):
        self.data = vel
        self.vel  = vel
        self.dt   = dt
        self.m    = m
        return

    def _window(self, W):
        n = self.vel.shape[0]
        w = [W(i,n) for i in range(0,n)]
        
        self.window = np.reshape(w, (n,1,1))
        self.data  *= self.window
        return

    def _padZeros(self, f):
        i = self.vel.shape[0]
        j = self.vel.shape[1]
        k = self.vel.shape[2]
        z = np.zeros((i*f,j,k))
        
        z[:i]       = self.data
        self.data   = z
        return

    def _vdos(self):
        k = 29979245800.0 # Divide by this to convert from Hz to cm^-1
        I = fft(self.c)
        v = fftfreq(len(self.c), self.dt) / k
        I = I[:len(v)//2]
        v = v[:len(v)//2]

        self.v = np.abs(v)
        self.I = np.abs(I)
        return

    def getSpectrum(self, win, pad, mir):
        if win:
            self._window(win)
        
        if pad:
            self._padZeros(pad)

        if mir:
            self.data = np.hstack((np.flipud(self.data), self.data))

        self.c = vacf(self.data, self.m)
        
        self._vdos()
        return
