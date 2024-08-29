import sys
import os
import numpy as np
import glob

from AaronTools.finders import BondedElements
from AaronTools.internal_coordinates import Bond
from AaronTools.geometry import Geometry
from AaronTools.fileIO import FileReader, read_types
from AaronTools.utils.utils import get_filename, glob_files
from AaronTools.theory import *

theory = Theory(
    method="MN15L",
    basis="def2-SVP",
    processors=12,
    memory=24,
    route={"scf": ["xqc", "MaxCycles=500"], "opt": ["MaxCycles=500"], "Freq": ["NoRaman"]}
)

for f in glob.glob("*.log", recursive=True):
    infile = FileReader(f, just_geom=False)
    geom = Geometry(infile, refresh_connected=True, refresh_ranks=False)

    theory.charge = infile["charge"]
    theory.multiplicity = infile["multiplicity"]
    geom.write(outfile=os.path.join(f"{f.split('.')[0]}.com"), theory=theory)


