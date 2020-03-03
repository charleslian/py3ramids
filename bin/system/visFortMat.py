#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Aug 11 21:29:16 2018

@author: clian
"""
import matplotlib
import numpy as np
from matplotlib import pyplot as plt
from py3ramids.plot.PlotUtility import scanFolder
import py3ramids.plot.setting as ma
import argparse

def read(filename):
    lines = open(filename).readlines()
    frames = [(i, line) for i, line in enumerate(lines) if '----' in line]
    data = []
    testLine = frames[0][0]+1
    if '(' in lines[testLine]:
        dataType = complex
    else:
        dataType = float
    for i, line in frames:
        ncol, nrow = [int(i) for i in line.split()[-2:]]
        nline = ncol*nrow
        dataFrame = np.array([dataType(line.replace(' ','').replace('+-','-')) for line in lines[i+1:i+1+nline]])
        data.append(dataFrame.reshape([nrow,ncol]))
    data = np.array(data)
    return data, dataType

# start the plotting code 
if __name__ == '__main__':
    class options(argparse.ArgumentParser):
        def __init__(self):
          super(options, self).__init__()
          
          self.add_argument('file', type=str, nargs='+',
                            help='the matrix data to visulize')
          self.add_argument('-w', '--write', type=bool, nargs='?',
                            help='write the data in numpy format')                
          self.add_argument('-p', '--plot', type=bool, nargs='?',
                            help='plot the data')      
          self.args = self.parse_args()
    
    options = options()
    args = options.args
    
    
    #print(args.write)
    
    for filename in args.file:
        print('Processing %s'%filename)
        if args.write:
            data, dataType = read(filename)
            np.save(filename,data)
        fig, axs = plt.subplots(data.shape[0],1,sharex=True,sharey=False,figsize=(6,data.shape[0]*3.5))
        if dataType is complex:
            vData = np.absolute(data)
        else:
            vData = data
        for i, frame in enumerate(vData):
            if data.shape[0] ==1:
                ax = axs
            else:
                ax = axs[i]
            vmax = frame.max()
            nRow, nCol = frame.shape
            ct = ax.imshow(frame,cmap='binary',vmax=vmax,vmin=0,aspect='auto',extent=[1, nCol+1, 0.5, nRow+0.5])
            plt.colorbar(ct, ax=ax)
            kargs=ma.getPropertyFromPosition(xlabel='',
                                             ylabel=r'',
                                             title='')
            ma.setProperty(ax,**kargs)
        
        fig.align_ylabels(axs)
        plt.tight_layout()
        
        #plt.axis('equal')
        SaveName = filename
        if args.plot:
          for save_type in ['.png']:
            filename = SaveName + save_type
            plt.savefig(filename,dpi=600,transparent=True)
