#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun May 14 14:57:18 2017

@author: clian
"""
import pyramids.io.result as dp
import pyramids.plot.setting as ma
import pyramids.plot.PlotUtility as ppu
from pyramids.io.fdf import tdapOptions
from matplotlib import pyplot as plt
import numpy as np
import os

Ry = 13.6
f = os.popen('grep "!    total energy" result')
energy = np.array([float(line.split()[-2]) for line in f.readlines()])
print energy.shape
fig, axs = plt.subplots(3,1,sharex=True,figsize=(6,8))
axs[0].plot((energy[3:]-energy[0])*Ry)

Efield = [[float(i) for i in line.split()] for line in open('TDAFIELD')]
Efield = np.array(Efield)/1E5
print Efield.shape
axs[1].plot(Efield)



def readData(filename='pwscf.phase.dat'):
  f = open(filename)
  text = f.readlines()
  nbnd, nkstot = [int(i) for i in text[0].split()]
  kweight = [float(i) for i in text[1].split()]
  
  nstep = (len(text) - 2)/(nkstot+1)
  #print nbnd, nkstot, nstep
  del text[1]
  del text[::(nkstot+1)]
  #print text
  #data = np.zeros([nbnd,nkstot,nstep])
  #data = 
  data = np.array([[float(i) for i in line.split()] for line in text])
  data = data[:nstep*nkstot].reshape([nstep, nkstot, nbnd])
  #print data[0,0,:]
  #[ for k in range(nkstot) for step in range(nstep)]
  return data, np.array(kweight)


#time, msd = dp.readMSD()
nocc = int(float(os.popen('grep "number of electrons" result').readline().split()[-1])/2.0)
print nocc

norm, kweight = readData('pwscf.norm.dat')
excite = (norm - norm[0,:,:])
for ib in range(excite.shape[2]):
    excite[:,:,ib] *= kweight
    #pass
time = np.array(range(excite.shape[0]+1))
step = min(time.shape[0]-1, excite.shape[0])
axs[2].plot(time[1:step+1], (excite[:step,:,:nocc]).sum(axis=(1,2)))
axs[2].plot(time[1:step+1], (excite[:step,:,nocc:]).sum(axis=(1,2)))