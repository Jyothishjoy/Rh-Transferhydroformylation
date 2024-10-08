#!/usr/bin/env python3

import sys
import os
import numpy as np
import glob

from AaronTools.finders import BondedElements,NotAny,BondedTo
from AaronTools.internal_coordinates import Bond
from AaronTools.geometry import Geometry
from AaronTools.fileIO import FileReader, read_types
from AaronTools.utils.utils import get_filename, glob_files
from AaronTools.component import Component
from AaronTools.theory import *

theory_I = Theory(
    method="M06L",
    basis="def2-TZVP",
    processors=12,
    memory=24,
    route={"scf": ["xqc", "MaxCycles=100"], "POP": ["NBO,AllOrbitals,ThreshOrbitals=5"], "NMR":["GIAO"]}
    )

theory_II = Theory(
    method="M06L",
    basis="def2-TZVP",
    processors=12,
    memory=24,
    route={"scf": ["xqc", "MaxCycles=100"], "pop": ["NBO,AllOrbitals,ThreshOrbitals=5"]}
    )

for f in glob.glob("*.log", recursive=True):
    dirname = os.path.dirname(f)
    infile = FileReader(f, just_geom=False)
    geom = Geometry(infile, refresh_connected=True, refresh_ranks=False)
    metal = geom.find("Pd")[0]
    P1 = geom.find("P")[0]
    P2 = geom.find("P")[1]
    
    ligand = geom.get_fragment(P1, stop=metal, as_object=True)

    theory_I.charge = 0
    theory_I.multiplicity = 1
    geom.write(outfile=os.path.join(dirname, f"{f[:-4]}_NBO.com"), theory=theory_I)
    ligand.write(outfile=os.path.join(dirname, f"{f[:-4]}_Lig_NBO.com"), theory=theory_I)
    
    theory_II.charge = 1
    theory_II.multiplicity = 2
    ligand.write(outfile=os.path.join(dirname, f"{f[:-4]}_Lig_Cat_NBO.com"), theory=theory_II)
    
    theory_II.charge = -1
    theory_II.multiplicity = 2
    ligand.write(outfile=os.path.join(dirname, f"{f[:-4]}_Lig_Ani_NBO.com"), theory=theory_II)
    
    
    