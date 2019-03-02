#! /usr/bin/env python3

from ase.io.trajectory import Trajectory
traj = Trajectory('qn.traj')
atoms= traj[-1]

print(atoms)
print('positions=',atoms.positions)


