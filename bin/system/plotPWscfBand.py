#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon May 29 11:06:46 2017

@author: clian
"""
import numpy as np
import argparse
  
class options(argparse.ArgumentParser):
    def __init__(self):
        super(options, self).__init__()
        self.add_argument('-yl', '--ylimits', type=float, nargs=2, default=None,
                          help='lower and upper limit of y')
        self.add_argument('-ef', '--fermiEnergy', type=float, nargs=1, default=0.0,
                          help='Fermi Energy')
        self.args = self.parse_args() 
    
options = options()
args = options.args
print(args)
#def readData(): 
#    f=np.loadtxt('bands.out.gnu')
#    data = []
#    nband = 0
#    for line in f.readlines():
#        if len(line) != 1:
#            data.append([float(i) for i in line.split()])
#        else:
#            nband += 1
#    data = np.array(data)
#    nkpt = data.shape[0]//nband 
#    band = np.zeros([nkpt,nband+1])
#    band[:,0] = data[:nkpt,0]
#    band[:,1:] = data[:,1].reshape([nband,nkpt]).T
#    #print band
#    return band

data = np.loadtxt('bands.out.gnu')#readData()
#efermi = 0
#data[:,1:] -= efermi

import matplotlib.pyplot as plt
import pyramids.plot.setting as ma
SaveName = __file__.split('/')[-1].split('.')[0]

fig, axs = plt.subplots(1,1,sharex=False, sharey=True,figsize=(8,6))#,figsize=(10,6
ax = axs
ax.plot(data[:,0],data[:,1:]-args.fermiEnergy,'.b', mfc='w',ms=3.0)
args = ma.getPropertyFromPosition(#0, 
                                  #title = 'Band',
                                  ylimits= args.ylimits,
                                  ylabel = 'Energy (eV)', 
                                  #xlimits=[0.9,1.2], 
                                  xlimits=[data[:,0].min(),data[:,0].max()], 
                                  #ylimits=[-0.7,0.4], 
                                  #ylimits=[-15, 5], 
                                  #xticklabels=[],
                                  #hline=[0.0],
                                  hline=[0],
                                  )
ma.setProperty(ax,**args)
plt.tight_layout()


plt.tight_layout()
for save_type in ['.pdf']:
  filename = SaveName + save_type
  plt.savefig(filename,transparent=True,orientation='portrait',dpi=600)


