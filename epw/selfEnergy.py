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
import os

def getData(filename = 'epw.pristine.out', process=0):
    if not os.path.exists('%s.npy'%filename) or True:
        output = open(filename).readlines()
        seLines = [i for i, line in enumerate(output) if 'Re[Sigma]' in line]
        memLines = [i for i, line in enumerate(output) if 'Memory usage:' in line][-1]
        kLines = [line for i, line in enumerate(output) if 'ik =' in line]
        kcoors = np.array([[float(j) for j in line.split()[-3:]] for line in kLines])
        startLine = seLines[-1]
        endLine = memLines
        data2=np.loadtxt(output[startLine+3:endLine-2])
        np.save(filename,data2)
    else:
        data2 = np.load('%s.npy'%filename)
    if process == 1:
      bands = set(data2[:,1])
      nbnd = len(bands)
      data = data2.reshape([nbnd, data2.shape[0]//nbnd, data2.shape[1]])
    elif process == 2:
      seq = np.argsort(data2[:,2])
      data = data2[seq,:]
    else:
      data = data2
    return data
def plot(axs, filename = 'epw.pristine.out'):
    titles = 'ik', 'ibnd', 'enk (eV)', 'Re[Sigma] (meV)', 'Im[Sigma] (meV)',  'Z', 'lam'
    
    dataAll = getData(filename, process=2)
    
    ax = axs
    ax.plot(dataAll[:,2], dataAll[:,4],'o-',lw=3,label=filename)

    kargs=ma.getPropertyFromPosition(ylabel=r'Im($\epsilon$) (meV)', xlabel='Energy (eV)',
                                     ylimits=[0,None], 
                                     xlimits=[-2,2],
                                     vline=[0])
    ma.setProperty(axs,**kargs)
    
def plotLineWidth(ax, filename, modes=None, **kargs):
    data = np.loadtxt(filename)
    all_modes = set(data[:,3])
    nmodes = len(all_modes)
    all_bands = set(data[:,1])
    nbnd = len(all_bands)
    #print(data.shape)
    #print(data[:,1:5].reshape([-1,nmodes,4]))
    data_ref = data[:,1:5].reshape([-1,nmodes,4])
    #print() 
    #print(data)
    for i in range(nmodes):
      if modes != None and i not in modes:
        continue
      else:
        mode = data_ref[:,i,:]
        seq = np.argsort(mode[:,1])
       
        ax.plot(mode[seq,1], mode[seq,3],'-',label=str(i+1))

    kargs=ma.getPropertyFromPosition(ylabel=r'Im($\epsilon$) (meV)', xlabel='Energy (eV)',
                                     ylimits=[0,None], 
                                     xlimits=[-2,2], #title=filename,
                                     vline=[0])
    ma.setProperty(ax,**kargs)
    
def getPrtgkk(filename, restart = False, **kargs):
    saveName = filename+'.prtg'
    if not os.path.exists('%s.npy'%saveName) or restart:
        output = open(filename).readlines()
        qLines = [line for i, line in enumerate(output) if 'iq =' in line]
        kLines = [line for i, line in enumerate(output) if 'ik =' in line]
        #kcoors = np.array([[float(j) for j in line.split()[-3:]] for line in kLines])
        #startLine = seLines[-1]
        #endLine = memLines
        #data2=np.loadtxt(output[startLine+3:endLine-2])
        #np.save(filename,data2)
    else:
        data2 = np.load('%s.npy'%saveName)
      
    print(data2)
    
def readBand(filename):
    text = open(filename).readlines()
    band = np.loadtxt(text[2::2])
    cood = np.loadtxt(text[1::2])
    print(cood)
    x = np.zeros(band.shape[0])
    x[0] = 0
    x[1:] = np.array([np.linalg.norm(cood[i] - cood[i-1]) for i in range(1, x.shape[0])])
    x = np.add.accumulate(x)
    #print(x)
    return x, band