import pytest, sys
import numpy as np

def test_coulomb():
    sys.path.append('../sourcecode/')
    from MvH_CO_JM8 import V_coul

    E, F       = V_coul([1,0,0], 1)
    assert  E == 1
    assert (F == [1,0,0]).all()

    E, F       = V_coul([1,1,0], 1)
    assert  E == 1 /  np.sqrt(2)
    x          = 1 / (np.sqrt(2))**3
    print(F)
    assert (F == [x,x,0]).all()

    E, F       = V_coul((1,1,1), 1)
    assert  E == 1 /  np.sqrt(3)
    x          = 1 / (np.sqrt(3))**3
    assert (F == [x,x,x]).all()

    return
