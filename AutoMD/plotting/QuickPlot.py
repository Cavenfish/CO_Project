from ..config import *

class QuickPlot:
    '''
    Quick plots for energy dissipation study
    '''
    #Default plot variables
    labels = {'co':'CO', '13co':r'$^{13}$CO', 'c18o':r'C$^{18}$O',
              '13c18o':r'$^{13}$C$^{18}$O'}

    colors = {'co':(1,0,1,1), '13co':(0,0.75,0.75,1),
              'c18o':(0.94,0.31,0.15,1), '13c18o':(0.27,0.04,0.16,1)}

    def __init__(self, h5):
        self.hdf = pd.HDFStore(h5)

    def _labelPlot(self):
        plt.xlabel('Time (ps)',   fontsize=16)
        plt.ylabel('Energy (eV)', fontsize=16)
        plt.legend()
        plt.tight_layout()

    def _leaves(self, s):
        for i in range(2):
            node = self.hdf.get_node(s + f'{i}')
            for j in node.__members__:
                if 'tau' in j:
                    continue
                yield f'{s}{i}/{j}'

    def _prep(self, df, col, n):
        x = df['Time']

        if n == 0:
            y = df[col]
        else:
            y = savgol_filter(df[col], n, 2)
        return x,y

    def spaget(self, dName, col='Sliced Energy', n=51):
        colors = [(i,0,i) for i in np.arange(1,0.15,-0.00425)]
        i=0
        for leaf in self._leaves(dName):
            df  = self.hdf.get(leaf)
            x,y = self._prep(df, col, n)
            plt.plot(x, y, color=colors[i])
            i  +=1

        df  = self.hdf.get(dName + 'Avg')
        x,y = self._prep(df, col, n)
        plt.plot(x, y, '--', color='gold', label='Average', lw=1)

        self._labelPlot()
        return plt.gca()

    def isotopes(self, dName, cluster, col='Sliced Energy', n=0):
        isos = ['co', '13co', 'c18o', '13c18o']
        for iso in isos:
            df  = self.hdf.get(dName + f'{iso}/{cluster}Avg')
            x,y = self._prep(df, col, n)
            plt.plot(x, y, label=self.labels[iso], color=self.colors[iso])

        plt.ylim((-0.05,0.45))
        plt.xlim((0,2000))
        plt.yticks([0,0.1,0.2,0.3,0.4])
        plt.xticks([0,500,1000,1500,2000])
        self._labelPlot()
        return plt.gca()

    def compare(self, names, labels, colors, col='Sliced Energy', n=51):
        for i in range(len(names)):
            df  = self.hdf.get(names[i])
            x,y = self._prep(df, col, n)
            plt.plot(x, y, label=labels[i], colors=colors[i])

        self._labelPlot()
        return plt.gca()

    def vacfOneTime(self, dName, time):
        c = {'eI': '#8EE3EF', 'r1': '#9055A2', 'r2': '#D499B9', 'r3': '#EACBD2'}

        df = self.hdf.get(dName)
        v  = df.index.values
        fI = df[f'{time} Full'    ].values
        eI = df[f'{time} (0, 1)'  ].values
        r1 = df[f'{time} (1, 10)' ].values
        r2 = df[f'{time} (10, 20)'].values
        r3 = df[f'{time} (20, 30)'].values

        plt.plot(v, fI,       label='Full',    color='black')
        plt.plot(v, eI, '--', label='Excited', color=c['eI'])
        plt.plot(v, r1, '--', label='r1',      color=c['r1'])
        plt.plot(v, r2, '--', label='r2',      color=c['r2'])

        if not np.array_equal(fI, r3):
            plt.plot(v, r3, '--', label='r3', color=c['r3'])

        return plt.gca()
