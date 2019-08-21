#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 24 08:59:33 2018

@author: clian
"""

from ase.dft.kpoints import get_special_points
from ase.build import bulk
from ase.io import read


#si = bulk('Si', 'diamond', a=5.459)

atoms = read('result',format='espresso-out')
points = get_special_points(atoms.cell)


line = 'K_POINTS crystal_b\n'
line += '%i\n'%len(points)
denK = 40
for point in points:
    line += '%4.3f %4.3f %4.3f %i ! %s \n' %(points[point][0],points[point][1],points[point][2],denK, point) 
    #print(points)


open('bandkpt.part','w').write(line)