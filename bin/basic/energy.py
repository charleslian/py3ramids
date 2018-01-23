#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec  5 21:04:25 2017

@author: clian
"""



def plot(ax, var):
  ax.plot(var.time, var.energy[0], '-',label='Total')
  ax.plot(var.time, var.energy[1], '--', label='KS')
  import py3ramids.plot.setting as ma
  ma.setProperty(ax, ylabel='Energy (eV)', xlabel='Time (fs)')


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
  fig, ax = plt.subplots(1,1)
  
  plot(ax, var)
  
  plt.tight_layout()
  SaveName = __file__.split('/')[-1].split('.')[0]
  plt.savefig(SaveName+'.pdf')
  
  

