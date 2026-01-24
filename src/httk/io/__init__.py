from httk.core import load, register_loader
from .cif.cif_reader import read_cif  # importing registers the loader

__all__ = ["read_cif"]

