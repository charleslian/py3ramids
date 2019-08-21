#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Aug 11 21:29:16 2018

@author: clian
"""

import numpy as np
from matplotlib import pyplot as plt
from py3ramids.plot.PlotUtility import scanFolder
import py3ramids.plot.setting as ma


# start the plotting code 
def plot(iv, ax):
    label = ''
    ax.plot([],[],label=label)
    kargs=ma.getPropertyFromPosition(xlabel='',
                                     ylabel=r'',
                                     title='')
      
      
    ma.setProperty(ax,**kargs)
    ma.add_label(iv, ax)
    columes = (4,5,6,10)
    lineK = np.array([[float(v) for i, v in enumerate(line.replace(')',' ').split()) if i in columes] 
            for line in open('bandK')])
    
    lineK[:,-1] = 1E-6
    meshK = np.array([[float(v) for i, v in enumerate(line.replace(')',' ').split()) if i in columes] 
            for line in open('meshK')])
    
    lines = 'K_POINTS tpiba\n'
    
    print(meshK.shape[0], lineK.shape[0])
    lines += '%i\n'%(meshK.shape[0] + lineK.shape[0])
    
    for points in [meshK, lineK]:
        for point in points:
            print(point)
            lines += '%12.6f\t%12.6f\t%12.6f\t%12.6f\n'%(point[0],point[1],point[2],point[3])
    print(lines)
    open('kout','w').write(lines)
    #data[:,-1] = 1E-5
    
    
    #print(data)
    #print(np.loadtxt(file, usecols=(4,5,6)))
    return kargs
# end the plotting code 
    


if __name__ == '__main__':
    fig, axs = plt.subplots(1,1,sharex=False,sharey=False,figsize=(6,8))
    plot(0,axs)
    plt.tight_layout()
    fig.align_ylabels(axs)
    SaveName = __file__.split('/')[-1].split('.')[0]
    if True:
      for save_type in ['.pdf','.png']:
        filename = SaveName + save_type
        plt.savefig(filename,dpi=600)