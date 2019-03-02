#!/usr/bin/env python3

from ase.io import write, read

# Part 2: Spectra calculations
#import ase.build
#ase.build.niggli_reduce(atoms)
from ase.io.trajectory import Trajectory
traj = Trajectory('qn.traj')
atoms = traj[-1]

write('POSCAR.vasp',atoms)
