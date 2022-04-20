import pytest, sys
sys.path.append('../')
sys.path.append('./')
import AutoMD as md

def test_vib_excitation():
    #10-step NVT at start gives initial momenta array

    xyz = 'tmp.xyz'
    xyz = md.run_langevinMD(xyz, n=10)
    xyz = md.excite_molecule(xyz, [0,1], 1.7)
    sys = md.read(xyz)
    pos = sys.get_positions()[:2]
    mas = sys.get_masses()[:2]
    vel = sys.get_velocities()[:2]
    E   = md.calc_Evib(pos, mas, vel)
    print(E)
    assert (E >= 1.6 and E <= 1.8)
    return

def test_rot_excitation():
    xyz = 'tmp.xyz'
    xyz = md.run_langevinMD(xyz, n=10)
    xyz = md.rotat_excite(xyz, [0,1], 0.0002)
    sys = md.read(xyz)
    pos = sys.get_positions()[:2]
    mas = sys.get_masses()[:2]
    vel = sys.get_velocities()[:2]
    E   = md.calc_Erot(pos, mas, vel)
    print(E)
    assert (E >= 0.0001 and E <= 0.0003)
    return

def test_trans_excitation():
    xyz = 'tmp.xyz'
    xyz = md.run_langevinMD(xyz, n=10)
    xyz = md.trans_excite(xyz, [0,1], 1)
    sys = md.read(xyz)
    pos = sys.get_positions()[:2]
    mas = sys.get_masses()[:2]
    vel = sys.get_velocities()[:2]
    E   = md.calc_Etran(pos, mas, vel)
    print(E)
    assert (E >= 0.9 and E <= 1.1)
    return
