from ..config import *

#Single exponential decay
def f(x, E, tau):
    return (E * np.exp(-x / tau))

#Two channel exponential decay
def g(x, E, tau, E2, tau2):
    return f(x, E, tau) + f(x, E2, tau2)

#The two channel expoential decay fits the
#decay curves best. The assumption here is:
#-> One channel is for vib-vib transfer
#-> Other is for vib-lattice transfer

def leaves(s):
    node = hdf.get_node(s + f'{i}')
    for j in node.__members__:
        yield f'{s}{i}/{j}'

def half_life(h5):
    hdf = pd.HDFStore(h5)

    for i in hdf.walk():
        tmp = {'uid':[], 'E1':[], 'tau1':[], 'E2':[], 'tau2':[]}
        for j in i[-1]:
            df = hdf.get(f'{i[0]}/{j}')
            p0 = (0.3, 100, 0.1, 10)

            popt, pcov = curve_fit(df['Time'], df['Sliced Energy'], p0)

            tmp['uid' ].append(j)
            tmp['E1'  ].append(popt[0])
            tmp['tau1'].append(popt[1])
            tmp['E2'  ].append(popt[2])
            tmp['tau2'].append(popt[3])

        hdf.put(f'{i[0]}/tau', tmp)
    return
