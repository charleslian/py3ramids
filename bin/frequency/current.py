#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec  5 21:04:25 2017

@author: clian
"""



def plot(ax, var):
  directions = ['x','y','z'] 
  linestyle = '-','--','-.'
  from scipy.fftpack import fft
  import numpy as np
  for idir in range(3):
    length = var.current.shape[0]//2
    #print(length)
    beta = 1E-3
    fftc = fft(var.current[:,idir]*np.exp(-beta*var.time))
    resE = 4.1356/var.time[-1]
    energy = np.arange(length)*resE #/ 2
    
    ax[1].plot(energy, np.abs(np.real(fftc[:length])), ls=linestyle[idir], label=directions[idir])
    ax[0].plot(energy, -np.imag(fftc[:length]), ls=linestyle[idir], label=directions[idir])
  import py3ramids.plot.setting as ma 
  ma.setProperty(ax[1], ylabel='Im $\epsilon$ (a.u.)', xlabel='Energy (eV)',
                 ylimits=[0,0.02], #xlimits=[0,8]
                 )
  ma.setProperty(ax[0], ylabel='Re $\epsilon$ (a.u.)', xlabel='Energy (eV)',
                 ylimits=[-0.02,0.02], 
                 #xlimits=[0,10]
                 )
  
if __name__ == '__main__':
  import os
  if os.path.exists('input.fdf'):
    soft = 'tdap'
  elif os.path.exists('input.in'):
    soft = 'tdpw'
  
  if soft == 'tdap':
    from py3ramids.io.SIESTA_interface import TdapVarible
    var = TdapVarible()
  else:
    from py3ramids.io.ESPRESSO_interface import TdpwVarible
    var = TdpwVarible()
  
  import matplotlib.pyplot as plt
  fig, ax = plt.subplots(2,1,figsize=(8,8),sharex=True)
  
  plot(ax, var)
  
  plt.tight_layout()
  SaveName = __file__.split('/')[-1].split('.')[0]
  plt.savefig(SaveName+'.pdf')
  
  

