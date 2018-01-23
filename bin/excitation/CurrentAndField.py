#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat May 13 07:43:17 2017

@author: clian
"""

import pyramids.io.result as dp
import pyramids.plot.setting as ma
import pyramids.plot.PlotUtility as ppu
from pyramids.io.fdf import tdapOptions
from matplotlib import pyplot as plt
import numpy as np
import os

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

fig, axs = plt.subplots(2,2,sharex=True,figsize=(8,5))
axs = axs.flatten()
norm, kweight = readData('pwscf.norm.dat')
excite = (norm - norm[0,:,:])
for ib in range(excite.shape[2]):
    excite[:,:,ib] *= kweight
    #pass
time = np.array(range(excite.shape[0]+1))
step = min(time.shape[0]-1, excite.shape[0])
axs[0].plot(time[1:step+1], (excite[:step,:,:nocc]).sum(axis=(1,2)))
axs[0].plot(time[1:step+1], (excite[:step,:,nocc:]).sum(axis=(1,2)))
args = ma.getPropertyFromPosition(0,title='Excited electrons',ylabel='n (e)')
ma.setProperty(axs[0],**args)
Efield = [[float(i) for i in line.split()] for line in open('TDEFIELD')]
Efield = np.array(Efield)/1E5

ppu.plotTotalEnergy(axs[2])
ppu.plotEField(axs[1])


f = os.popen('grep "current is " result')
current = np.array([[float(i) for i in line.split()[-3:]] for line in f.readlines()])
for idir in range(3):  
    axs[3].plot(time[1:], current[:,idir],'-',label=['x','y','z'][idir])

args = ma.getPropertyFromPosition(3,title='Current',ylabel='j (a.u.)')
ma.setProperty(axs[3],**args)

plt.tight_layout()

SaveName = __file__.split('/')[-1].split('.')[0]
for save_type in ['.png']:
  filename = SaveName + save_type
  plt.savefig(filename,orientation='portrait',dpi=600)