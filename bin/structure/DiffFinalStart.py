#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  7 08:55:46 2017

@author: cl-iop
"""

import os
from ase.visualize import view
from ase.io.trajectory import Trajectory
traj = Trajectory('Trajectory')
import py3ramids.plot.PlotUtility as ppu


start, final = traj[1], traj[-1]
posDiff = final.positions - start.positions 


#start
demoImag = start.copy()
demoImag.set_velocities(posDiff)

ppu.generateStructPNG(start,filename='start',showcell=2, repeat = [2,2,1])
#ppu.generateStructPNG(start)
#print(posDiff)
#view(demoImag)