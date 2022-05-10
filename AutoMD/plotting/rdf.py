from ..config import *
from ase.geometry.analysis import Analysis

def rdf(xyz, rmax, nbins, elem):
    system = read(xyz)
    system.set_cell([100,100,100])
    
    ana = Analysis(system)
    g,r = ana.get_rdf(rmax=rmax, nbins=nbins, elements=elem, 
                      return_dists=True)[0]
    plt.plot(r, g)
    plt.show()
