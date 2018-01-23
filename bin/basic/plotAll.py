#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec  5 21:04:25 2017

@author: clian
"""



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
    
  
  
  #print(basic)
  import matplotlib.pyplot as plt
  fig, axs = plt.subplots(5,1,figsize=(6,12),sharex=True)
  
  import py3ramids.bin.basic.current
  py3ramids.bin.basic.current.plot(axs[0],var)
  
  import py3ramids.bin.basic.energy
  py3ramids.bin.basic.energy.plot(axs[1],var)
  
  import py3ramids.bin.basic.Afield
  py3ramids.bin.basic.Afield.plot(axs[2],var)
  
  import py3ramids.bin.basic.Efield
  py3ramids.bin.basic.Efield.plot(axs[3],var)
  
  import py3ramids.bin.basic.carrier
  py3ramids.bin.basic.carrier.plot(axs[4],var) 
  plt.tight_layout()
  SaveName = __file__.split('/')[-1].split('.')[0]
  plt.savefig(SaveName+'.pdf')
  
  

