import pytest
import sys

def test_numpy_import():
    import numpy
    return

def test_ase_import():
    import ase
    return

def test_MvH_import():
    sys.path.append('../sourcecode/')
    from MvH_CO_JM8 import MvH_CO
    return

