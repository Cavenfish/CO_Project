import pytest, sys
sys.path.append('./AutoMD/')
sys.path.append('../AutoMD/')
sys.path.append('../')
sys.path.append('./')

def test_numpy_import():
    import numpy
    return

def test_ase_import():
    import ase
    return

def test_MvH_import():
    from MvH_CO_JM8 import MvH_CO
    return

def test_package_import():
    import AutoMD as md
    return

