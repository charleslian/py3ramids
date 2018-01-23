#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 22 08:31:42 2017

@author: cl-iop
"""
import numpy as np
import pandas as pd
import os

def calculateDOS(var, step, xlimits=None, intepNum=1000, bins=200, ref=True, interp=True):
  
  kweight = var.kpts[:,3]
  x = np.arange(kweight.shape[0])
  
  xlabel = 'Energy (eV)'
  ylabel = 'Population'
  
  sortedDF = generatePopulationFile(var, step, ref)
  #print(sortedDF)
  xt = sortedDF[xlabel]
  yt = sortedDF[ylabel]
  
  #print 'yt is',yt  
  if xlimits != None: 
    x = xt[xt>xlimits[0]][xt<xlimits[1]]
    y = yt[xt>xlimits[0]][xt<xlimits[1]]
  else:
    x = xt
    y = yt
    
  dos, bin_edges = np.histogram(x, bins=bins, range=xlimits)
  par, bin_edges = np.histogram(x, bins=bins, range=xlimits, weights=y/2)
  #parDn, bin_edges = np.histogram(x,bins=bins,range=xlimits,weights=-y/2)
  if not interp:
    return bin_edges[:-1],dos,par
  
  def interp(xin,yin,xout):
    from scipy.interpolate import interp1d
    spline = interp1d(xin, yin, kind='cubic')
    return spline(xout) 
  
  eDosInterp = np.linspace(bin_edges[0], bin_edges[-2], intepNum)
  yDosInterp = interp(bin_edges[:-1], dos, eDosInterp)
  yParInterp = interp(bin_edges[:-1], par, eDosInterp)
  
  return eDosInterp, yDosInterp, yParInterp


def generatePopulationFile(var, step, ref=True):
  xlabel = 'Energy (eV)'
  ylabel = 'Population'
  

  #dataFilename = xlabel+'vs'+ylabel+'.csv'
  filename = xlabel+'vs'+ylabel+str(step)+'.csv'
  if os.path.exists(filename) and False:
    sortedDF =  pd.read_csv(filename)
  else:
    exe = var.readTDData(mode='norm')
    #print(exe.shape)
    eigen = var.readTDData(mode='value')
    #print(eigen)
    if ref:
      df = pd.DataFrame({xlabel:eigen[step].flatten(),
                         ylabel:(exe[step] - exe[0]).flatten()})
    else:
      df = pd.DataFrame({xlabel:eigen[step].flatten(),
                         ylabel:exe[step].flatten()})
    sortedDF = df.sort_values(by=xlabel)
    sortedDF.to_csv(filename)
  return sortedDF