"""
AutoMD package for doing data analysis
"""

from .half_life        import half_life
from .diffusion_length import diffusion_length
from .local_geo        import local_geo
from .unpickle         import unpickle
from .get_surface      import get_surface

__functions__ = [half_life, diffusion_length]
