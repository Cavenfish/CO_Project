from ..config   import *
from ..analysis import VACF

def slices_of_vacf(traj, frames, n):
    tj = Trajectory(traj)
    df = {}

    for i in frames:
        #prep system
        atms = tj[i]
        calc = MvH_CO(atoms=atms)
        atms.set_calculator(calc)
        dyn  = VelocityVerlet(atms, units.fs)
        tmp  = Trajectory('temp.traj', 'w', atms)
        dyn.attach(tmp.write, interval=3)

        #Run NVE
        dyn.run(n)

        #Compute VACF
        tmp = Trajectory('temp.traj')
        vel = np.array([i.get_velocities() for i in tmp])
        acf = VACF(vel, 3e-15)
        win = acf.HannMirror
        pad = 8
        mir = True
        acf.getSpectrum(win,pad,mir)

        #save VACF
        key     = f'{i/10} ps'
        df[key] = acf.I

    #Make DataFrame
    df = pd.DataFrame(df, index=acf.v)

    #Save dataframe
    df.to_csv('vacfTESTdata.csv')
    return df
