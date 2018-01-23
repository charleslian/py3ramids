from __future__ import print_function
import matplotlib.pyplot as plt
import py3ramids.plot.setting as ma


import ase.dft.band_structure as bs
from ase.io import read
atoms = read('POSCAR',format='vasp')

dftband = bs.BandStructure.read('band.json')
#data = dftband.todict()
from ase.dft.kpoints import bandpath
kpts = bandpath([[0.333,0.666,0.0],[0.5,0.0,0.0],[0.0,0.0,0.0],[0.666,0.333,0.0]],atoms.cell,npoints=90)


fig, ax = plt.subplots(1,1,sharex=True,figsize=(6,8))
x = kpts[-2]
y = dftband.energies[0] - 4.74503
ax.plot(x,y)
print(kpts[0])
ma.setProperty(ax, ylimits=[-2,2], vline=kpts[-1])
