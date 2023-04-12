import copy
from ..config import *
from numpy.fft import fftfreq, fft, ifft

class VACF:
    '''
    Velocity AutoCorrelation Function
    '''
    #Window Functions
    Hann       = lambda _,n,N: np.sin((np.pi*n)/N)**2
    Welch      = lambda _,n,N: 1 - ((n-N/2)/(N/2))**2
    HannMirror = lambda _,n,N: np.cos((np.pi*n)/(2*(N-1)))**2

    def  __init__(self, vel, dt):
        self.atms = False
        self.vel  = vel
        self.dt   = dt
        return

    def _window(self, W):
        n = self.c.shape[0]
        w = [W(i,n) for i in range(0,n)]

        self.window = w
        self.c     *= self.window
        return

    def _padZeros(self, f):
        i = len(self.c)
        z = np.zeros(i*f)

        z[:i]       = self.c
        self.c   = z
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

    def _vacf(self):
        T = self.vel.shape[0]
        N = self.vel.shape[1]
        D = self.vel.shape[2]
        c = np.zeros(T-1)

        if self.atms:
            atoms = self.atms
        else:
            atoms = range(N)
        
        for i in atoms:
            for j in range(D):
                data = self.vel[:,i,j]
                forw = fft(data, n=2*T)
                tmp  = forw * np.conjugate(forw)
                back = ifft(tmp)[:T] / T
                c   += np.real(back[:T-1])
 
        c   /= c[0]
        return c

    def getSpectrum(self, win, pad, mir):
        self.c = self._vacf()
        self.C = copy.deepcopy(self.c) 

        if win:
            self._window(win)

        if pad:
            self._padZeros(pad)

        if mir:
            self.c = np.hstack((np.flipud(self.c), self.c[1:]))

        self._vdos()
        return
