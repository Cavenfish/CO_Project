import pytest, sys
import numpy as np

def test_coulomb():
    sys.path.append('../AutoMD/')
    from MvH_CO_JM8 import V_coul

    E, F       = V_coul([1,0,0], 1)
    assert  E == 1
    assert (F == [1,0,0]).all()

    E, F       = V_coul([1,1,0], 1)
    assert  E == 1 /  np.sqrt(2)
    x          = 1 / (np.sqrt(2))**3
    assert (F == [x,x,0]).all()

    E, F       = V_coul((1,1,1), 1)
    assert  E == 1 /  np.sqrt(3)
    x          = 1 / (np.sqrt(3))**3
    assert (F == [x,x,x]).all()

    return

def test_exchange():
    sys.path.append('../AutoMD/')
    from MvH_CO_JM8 import V_exch
    import numpy as np

    #I think the radial vector used points in 
    #negative direction so the negative in the force
    #is dropped. Might want to look more into it.
    
    r_is_1 = np.array([1,0,0])

    E, F       = V_exch(r_is_1, 1, 0)
    assert  E == 1
    assert (F == [0,0,0]).all()

    E, F       = V_exch(r_is_1, 1, 1)
    assert  E == np.exp(-1)
    assert (F == [np.exp(-1),0,0]).all()

    #More comprehensive test with proper diff given
    c1 = np.array([1, 5,2])
    c2 = np.array([8,14,9])
    a  = b = 1

    diff = c2 - c1 
    
    #Note in sourcecode format we are getting F on c1
    #with radial vector going from c2 to c1
    
    E, F       = V_exch(diff, a, b) 
    assert  E == 1.5471621744977957e-6
    assert (F == 1.5471621744977957e-6 * diff/(np.sqrt(np.dot(diff,diff)))).all()

    return
