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
    return

def rdf_of_array(xyzs, rmax, nbins, elem):
    
    for xyz in xyzs:
        system = read(xyz)
        system.set_cell([100,100,100])

        ana = Analysis(system)
        g,r = ana.get_rdf(rmax=rmax, nbins=nbins, elements=elem,
                      return_dists=True)[0]
        plt.plot(r, g)
    plt.show()
    return

def rdf_of_traj(traj, rmax, nbins, elem):
    tmp  = Trajectory(traj)
    sys0 = tmp[0]
    sys1 = tmp[-1]
    sys0.set_cell([100,100,100])
    sys1.set_cell([100,100,100])

    ana0 = Analysis(sys0)
    ana1 = Analysis(sys1)

    g0,r0 = ana0.get_rdf(rmax=rmax, nbins=nbins, elements=elem,
                        return_dists=True)[0]
    g1,r1 = ana1.get_rdf(rmax=rmax, nbins=nbins, elements=elem,
                        return_dists=True)[0]

    plt.plot(r0,g0, color='blue')
    plt.plot(r1,g1, color='red')
    plt.show()
    return


