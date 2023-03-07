from ..config import *
from yaml import safe_load
from alphashape import alphashape
from shapely.geometry import MultiPoint

class Structure:
    def __init__(self, inp):
        with open(inp, 'r') as f:
            self.p = safe_load(f)

        self.system = read(self.p['xyz'])
        self.aShape = alphashape(self.system.positions)

    def getGlobalDensity(self):
        v = self.aShape.volume * 1e-24
        n = len(self.system.positions) // 2
        g = n * (28.0101 / 6.022e23)
        d = g / v

        units = 'g/cm3'
        print(f'Density: {d} {units}')
