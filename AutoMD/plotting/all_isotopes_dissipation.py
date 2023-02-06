from ..config import *
from operator import itemgetter

def all_isotopes_dissipation(isoDict, err=False, noBkg=False):
    def get_avg(csvNames):
        key = 'Sliced Energy'
        N   = len(csvNames)

        for fname in csvNames:
            df   = pd.read_csv(fname)
            
            try:
                avg += df[key].values
            except NameError:
                avg  = df[key].values
            except:
                N -= 1
                continue

        avg /= N
        time = df['Time']
        return avg, time

    def get_std(csvNames, avg):
        key = 'Sliced Energy'
        N    = len(csvNames)

        for fname in csvNames:
            df = pd.read_csv(fname)

            try:
                std += (df[key] - avg)**2
            except:
                std  = (df[key] - avg)**2
        
        std /= (N -1)
        std  = np.sqrt(std)
        return std

    # Define plot color variables
    labels   = {'co':'CO', '13co':r'$^{13}$CO',
                '13c18o':r'$^{13}$C$^{18}$O', 'c18o':r'C$^{18}$O'}
    
    colors   = {'co':(1,0,1,1), '13co':(0,0.75,0.75,1),
                '13c18o':(0.27,0.04,0.16,1), 'c18o':(0.94,0.31,0.15,1)}
    
    # Make plot
    for iso in isoDict:
        csvNames = isoDict[iso]
            
        avg, time = get_avg(csvNames)
        l,c       = labels[iso], colors[iso]

        if err:
            std = get_std(csvNames, avg)
            mks, cs, bs = plt.errorbar(time, avg, std, label=l, color=c,
                                       elinewidth=0.5, capsize=0.75)
            [bar.set_alpha(0.05) for bar in bs]
            [cap.set_alpha(0.05) for cap in cs]
        else:
            plt.plot(time, avg, label=l, color=c)


    plt.xlabel('Time (ps)',   fontsize=18)
    plt.ylabel('Energy (eV)', fontsize=18)
    #plt.title('Vibrational Energy Dissipation', fontsize=20)

    h, l = plt.gca().get_legend_handles_labels()
    tmp  = sorted(zip(h,l), key=itemgetter(1))
    h, l = zip(*tmp)
    plt.legend(h, l, fontsize=12)

    #return plot for any extra unique tweaks
    return plt.gcf()
