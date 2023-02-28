from ..config import *
from yaml import safe_load

class Structure:
    def __init__(self, inp):
        with open(inp, 'r') as f:
            p = safe_load(f)
        
        self.p = p
        self.system =
        
