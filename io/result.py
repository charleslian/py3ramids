#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 25 20:28:55 2016
The collection of the useful methods used in the TDAP results analysis
"""
__author__ = "Chao Lian <charleslian@126.com>"
__date__ = "Fri July 22 18:16:20 2016"

__version__ = "2.0"
__credits__ = "Chao Lian initial and maintain the codes"

import numpy as np
import os
from py3ramids.io.fdf import tdapOptions

def getElectronStepLength(filename='input.in'):
  for line in open(filename).readlines():
      if 'edt' in line.split():
          return float(line.split()[-1].replace(',',''))*0.04838

def getTrajactory():
  from ase.io.trajectory import Trajectory
  filename='Trajectory'
  return Trajectory(filename)
  
  
  
  
def getEELS(filename='q_list',prefix='EELS_'):
  Q = np.loadtxt(filename) 
  Z = []
  for i, q in enumerate(Q):
    filename = prefix + str(i+1) 
    d = np.loadtxt(filename, delimiter=',')
    E = d[:, 0]
    C = d[:, 2]
    Z.append(C)
  # to plot contour, transpose Z is needed (exchange x and y)
  contour = np.transpose(np.array(Z))
  X, Y = np.meshgrid(Q, E)
  
  return Q, E, X, Y, contour 

#------------------------------------------------------------------------------
def getTime():
  context = open('siesta.times','r').readlines()
  labels = ['Calls', 'Prg.com', 'Prg.com/Tot.com', 'Prg.tot', 'Prg.tot/Tot.tot', 'Nod.avg/Nod.max']
  #print labels

  data = [[float(num) for num in line.split()[1:]] for line in context[13:-13]]
  index = [line.split()[0] for line in context[13:-13]]
  #print index
  #print data
  import pandas as pd
  return pd.DataFrame(data = data, index = index, columns = labels)
#------------------------------------------------------------------------------

def getStructrue():
  if os.path.exists('structure.vasp'):
    from ase.io import read
    atoms = read('structure.vasp',format='vasp')
    #import ase.build
    #ase.build.niggli_reduce(atoms)
  elif os.path.exists('siesta.XV'):
    from ase.calculators.siesta.import_functions import xv_to_atoms
    atoms = xv_to_atoms('siesta.XV')
  return atoms

def pythonGrep(string,filename):
 import re
 return [line for line in open(filename, "r") if re.search(string, line)]

def recoverAllKPoints(kcoor, reciprocal_vectors, repeat = [1,1,1], flatten=True):
  #if plane == 'xy': 
  kall = []
  klist = []
  indexKList = np.arange(kcoor.shape[0])
  points = np.array([(reciprocal_vectors[0,:]*i+
                      reciprocal_vectors[1,:]*j+
                      reciprocal_vectors[2,:]*k) 
                      for i in range(-repeat[0],repeat[0]+1) 
                      for j in range(-repeat[1],repeat[1]+1)
                      for k in range(-repeat[2],repeat[2]+1)])
  for point in points:
    if flatten:  
      kall.extend(kcoor + point)
      kall.extend(-kcoor + point)
      klist.extend(indexKList)
      klist.extend(indexKList)
    else:
      kall.append(kcoor + point)
      kall.append(-kcoor + point)
      klist.append(indexKList)
      klist.append(indexKList)    
  return np.array(kall), np.array(klist)

  
def getImagElectricFunction(direction, dumping=0.00):
  from scipy.fftpack import fft
  timeArray, efield = getEField()  
  timeArray, dipole = getDipolePython()
  
  option = tdapOptions()
  lengthTime = option.tdTimeStep[0]
  
  
  timeArray = timeArray[2:]/option.tdTimeStep[0]
  numStep =  len(timeArray)
  
  dipole = dipole[2:,direction]
  dipole *= np.exp(-dumping*timeArray)
                            
  freqResolution = 1.0/(lengthTime*numStep)
  freqArray = (timeArray-(numStep/2.0))*freqResolution
  
  energyArray = freqArray*4.1356
  energyArray = energyArray[numStep/2:]
  
  epsilon = fft(dipole)[:numStep/2] 
  epsilon = (np.real(epsilon),np.imag(epsilon))
  import pandas as pd
  df = pd.DataFrame({'x-energy':energyArray, 
                    'y1-Im(alpha)':np.abs(epsilon[0]*energyArray/4.1356), 
                    'y2-Re(alpha)':np.abs(epsilon[1]*energyArray/4.1356)})
  #print df
  df.to_csv('DielectricFunction.csv',sep=',')
  return energyArray, np.abs(epsilon[0]*energyArray/4.1356), np.abs(epsilon[1]*energyArray/4.1356)
  
#-------------------------------------------------------------------
def loadSaved(filename):
  filename = filename + '.npy'
  if os.path.exists(filename):
    return np.load(filename)
  else:
    return []
#-------------------------------------------------------------------     
def readBerryFile(step = ''):
  """
  return the Eigenvalues read from systemLabel.EIG
  """
  options = tdapOptions()
  systemLabel = options.label
  BerryFile = open(systemLabel+step+'.Berry')
  Berry = [[float(value) for value in line.split()] for line in BerryFile.readlines()]
  return np.array(Berry)
#-------------------------------------------------------------------
def getBerry():
  options = tdapOptions()
  selectStep = getBerrySteps()
  #print selectStep
  timestep = options.tdTimeStep[0]
  selectTime = selectStep*timestep
  berry = np.array([readBerryFile(str(step)) for step in selectStep])
  return selectTime, berry
#
#-------------------------------------------------------------------  
def getBerrySteps():
  options = tdapOptions()
  stepLines = os.popen('ls '+options.label+'*.Berry').readlines()
  steps = [int(step.split('.')[0].replace(options.label,'')) for step in stepLines]
  return np.sort(steps)
  
def findBandPath(atoms, points, kLine, toleAngle = 0.001):
  reciprocal_vectors = 2*np.pi*atoms.get_reciprocal_cell()
  if os.path.exists('input.fdf'):
    kcoor, kweight = readKpoints()
    kall, klist = recoverAllKPoints(kcoor, reciprocal_vectors)
    #ax.plot(kcoor[:,0], kcoor[:,1],'o')
    #ax.plot(kall[:,0], kall[:,1],'.')
    
  selectedPoints = [(points[name], name) for name in kLine]
  
  lines = [(selectedPoints[i],selectedPoints[i+1]) 
            for i in range(len(selectedPoints) - 1)] 
  kpath = []     
  import pandas as pd
  
  last = 0.0
  xticks = []
  xticks.append(last)
  for a, b in lines:
    coorA = np.dot(a[0],reciprocal_vectors)
    coorB = np.dot(b[0],reciprocal_vectors)
    #ax.plot(coorB[0], coorB[1],'or',ms=20)
    #ax.plot(coorA[0], coorA[1],'or',ms=20)
    reVect =  coorB - coorA
    kpath_alongAB = []
    normReVect = np.linalg.norm(reVect) #AB
    for index, (k, ucellK) in enumerate(zip(kall, klist)):
      vect = k - coorA # AK
      normVect = np.linalg.norm(vect) 
      normProd = np.linalg.norm(vect)*np.linalg.norm(reVect)
      dotProd = np.dot(reVect, vect)
      #outProd = np.outer(reVect, vect)
      
      angleContained = (np.abs(dotProd - normProd) < toleAngle)
      if angleContained and normVect <= normReVect:
        #print np.linalg.det(outProd), ucellK
        kpath_alongAB.append((index, ucellK, dotProd/np.linalg.norm(reVect) + last))
    last += normReVect
    xticks.append(last)
    kpath.extend(kpath_alongAB)
   
  kpath = pd.DataFrame(kpath, columns=['sc-index','uc-index','distance'])
  kpath = kpath.sort_values(['distance'])
  #kpts = np.array(kpath.values[:,0],dtype=int)

  #kpath.to_csv('selectedKPath.csv')
  
  return xticks, kpath

  

#-------------------------------------------------------------------
def readEpsilonImaginary(filename = 'siesta.EPSIMG'):
  """
  return the imaginary of epsilon
  """
  epsimgFile = open(filename)
  epsimg = np.array([[float(num) for num in line.split()] 
  for il, line in enumerate(epsimgFile.readlines()) if il > 7])
  return epsimg
#-------------------------------------------------------------------     
def getFermiEnergy():
  """
  return the Fermi Energy read from systemLabel.EIG
  if failed, return 0 rather than raise an exception
  """
  options = tdapOptions()
  systemLabel = options.label
  if os.path.exists(systemLabel+'.EIG'):
    Efermi = float(open(systemLabel+'.EIG').readline().split()[0]) 
  else:
    Efermi = 0.0
  return Efermi
#-------------------------------------------------------------------

def getEnergyTemperaturePressure(ave=False):    
  """
  return the Temperature, KS Energy and Total Energy as the dimension of Nstep
  read from systemLabel.MDE
  returns the Temperature, KS Energy, Total Energy, Volume, Pressure
  """
  options = tdapOptions()
  systemLabel = options.label
  if options.mdTimeStep[0] < 1E-10:
    timestep = options.tdTimeStep[0]
    start = 2
  else:
    start = 0
    timestep = options.mdTimeStep[0]
  if(os.path.exists(systemLabel+'.MDE')):  
      energy_file = open(systemLabel+'.MDE')
      data = []
      for i in energy_file.readlines():
        if i.split()[0]=='#': continue
        data.append([float(j.replace('*','0')) for j in i.split()])
    
      data=np.array(data) 
      if ave:
        numAtom = getNumOfAtoms()
        #print numAtom
        data[start:,2:4] /= numAtom
      #print data
      X = data[:,0]*timestep
      return X[start:], data[start:,1], data[start:,2], data[start:,3], data[start:,4], data[start:,5]
  else:
      time, msd = readMSD()
      f = os.popen('grep "Ekin + Etot (const)" result')
      energy = np.array([float(line.split()[-2]) for line in f.readlines()])*13.6
      #data=np.array([energy,energy,energy,energy,energy]).T
      return time, energy, energy, energy, energy,energy
  
 #-------------------------------------------------------------------
def readMSD(filename='pwscf.msd.dat'):
  f = open(filename)
  text = f.readlines()
  data = np.array([[float(word) for word in line.split()] for line in text])
  return data[:,0]*1E3, data[:,1] 

#-------------------------------------------------------------------

#-------------------------------------------------------------------
def getAdiabaticEigenvalue(selectK=None):
  """
  """
  options = tdapOptions()
  selectStep = getEIGSteps()
  #print selectStep
  timestep = options.tdTimeStep[0]
  selectTime = selectStep*timestep
  
  SaveName = 'AdiabaticEigenvalue'
  exe = loadSaved(SaveName)
  if len(exe) != len(selectStep):
    exe = []
    for index,step in enumerate(selectStep):
      partition = readEigFile(options.label+str(step)+'.EIG')
      exe.append(partition)
    exe = np.array(exe)
    np.save(SaveName,exe)
  
  return selectTime, exe

#-------------------------------------------------------------------
def getProjectedPartition(selectK=None):
  """
  """
  options = tdapOptions()
  selectStep = getEIGSteps()
  #print selectStep
  timestep = options.tdTimeStep[0]
  selectTime = selectStep*timestep
  
  SaveName = 'ProjectedPartition'
  exe = []
  if len(exe) != len(selectStep):
    exe = []
    for index,step in enumerate(selectStep):
      partition = readEigFile(options.label+str(step)+'q.EIG')
      exe.append(partition)
    exe = np.array(exe)
    np.save(SaveName,exe)
  
  return selectTime, exe
#-------------------------------------------------------------------
def generatePopulationFile(step, ref=True):
  xlabel = 'Energy (eV)'
  ylabel = 'Population'
  
  import pandas as pd
  #dataFilename = xlabel+'vs'+ylabel+'.csv'
  filename = xlabel+'vs'+ylabel+str(step)+'.csv'
  if os.path.exists(filename) and False:
    sortedDF =  pd.read_csv(filename)
  else:
    time, exe = getProjectedPartition()
    time, eigen = getAdiabaticEigenvalue()
    if ref:
      df = pd.DataFrame({xlabel:eigen[step].flatten(),
                         ylabel:(exe[step] - exe[0]).flatten()})
    else:
      df = pd.DataFrame({xlabel:eigen[step].flatten(),
                         ylabel:exe[step].flatten()})
    
    
    sortedDF = df.sort_values(by=xlabel)
    sortedDF.to_csv(filename)
  return sortedDF
  
def calculateDOS(step, xlimits=None, intepNum=1000, bins=200, ref=True, interp= True):
  kcoor, kweight = readKpoints()
  x = np.arange(kcoor.shape[0])
  
  xlabel = 'Energy (eV)'
  ylabel = 'Population'
  
  sortedDF = generatePopulationFile(step,ref)
  
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


def getHomo():
  """
  return the index of the highest occupied molecular orbital (HOMO)
  """
  options = tdapOptions()
  #NumElectron=float(os.popen('grep "Total number of electrons:" result').readline().split()[-1]) 
  NumElectron=float(pythonGrep("Total number of electrons:",'result')[0].split()[-1]) 
  readFile = os.popen('grep -i "SpinPolarized" '+ options.inputFile)
  line = readFile.readline()
  if line == '':
    nspin = 2
  else:
    nspin = 1
  homo = int(NumElectron) / nspin
  return homo
#-------------------------------------------------------------------
def getDipoleFromFile(filename):
  lines = open(filename).readlines()
  dipole = [(float(line.split()[-3]),float(line.split()[-2]),float(line.split()[-1])) for line in lines]
  return np.array(dipole[:-1])
#-------------------------------------------------------------------
def getDipole():
  options = tdapOptions()
  lines = os.popen('grep -n "Electric dipole (a.u.)" result').readlines()    
  dipole = [(float(line.split()[-3]),float(line.split()[-2]),float(line.split()[-1])) for line in lines]
  timestep = options.tdTimeStep[0]
  dipole = np.array(dipole[:-1])
  time = np.arange(len(dipole))*timestep
  return time, dipole
  
#-------------------------------------------------------------------
def getDipolePython():
  options = tdapOptions()
  context = []
  for line in open('result').readlines():
      if len(line) == 1:
          continue
      #print line,
      if 'Electric' in line.split() and 'TDAP' in line.split() :
          context.append(line[:-1])
          
  #for line in context:
  #    print line.split()[-3:]
       
  dipole = [[float(line.split()[-3]),float(line.split()[-2]),float(line.split()[-1])] for line in context]
  #print context
  #dipole = [(float(line.split()[-3]),float(line.split()[-2]),float(line.split()[-1])) for line in lines]
  timestep = options.tdTimeStep[0]
  dipole = np.array(dipole)
  time = np.arange(len(dipole))*timestep
  return np.array(time), np.array(dipole)
#-------------------------------------------------------------------
def getIonicPolarization():
  lines = os.popen('grep -n "The ionic Polarization (a.u.):" result').readlines()    
  dipole = [(float(line.split()[-3]),float(line.split()[-2]),float(line.split()[-1])) for line in lines]
  return np.array(dipole) 
#-------------------------------------------------------------------
def getGaugeFreePolarization():
  lines = os.popen('grep -n "The gauge freedom Polarization" result').readlines()    
  dipole = [(float(line.split()[-3]),float(line.split()[-2]),float(line.split()[-1])) for line in lines]
  return np.array(dipole) 
#-------------------------------------------------------------------
def getElectronPolarization():
  lines = os.popen('grep -n "The electron Polarization (a.u.):" result').readlines()    
  dipole = [(float(line.split()[-3]),float(line.split()[-2]),float(line.split()[-1])) for line in lines]
  return np.array(dipole) 
#-------------------------------------------------------------------
def getEigAndPar(timeSteps=None, forceReread = False):
  """
  return the Times, Eigvalues and Partitions at selected timeSteps
  write two files EIG and PAR
  """
  options = tdapOptions()
  systemLabel = options.label
  timestep = options.tdTimeStep[0]

  if timeSteps == None:
    timeSteps = getEIGSteps()
  time = np.array(timeSteps)*timestep
  
  eig = loadSaved('EIG')
  par = loadSaved('PAR')
  saved = ((len(eig)) != 0 and (len(par)) != 0)
  if (not saved) or forceReread:
    eig = []; par = []
    for i in timeSteps:
      eig.append(readEigFile(systemLabel+str(i)+'.EIG'))
      par.append(readEigFile(systemLabel+str(i)+'q.EIG'))
    eig = np.array(eig)
    par = np.array(par)
    np.save('EIG',eig)
    np.save('PAR',par)
  
  return time, eig, par
#-------------------------------------------------------------------
def getTDEig(timeSteps=None, forceReread = False):
  """
  return the Times, Eigvalues and Partitions at selected timeSteps
  write two files EIG and PAR
  """
  options = tdapOptions()
  systemLabel = options.label
  timestep = options.tdTimeStep[0]
  steps = np.arange(3,options.mdFinalStep)
  
  time = np.array(steps)*timestep
  
  eig = loadSaved('TDEIG')
  saved = ((len(eig)) != 0)
  if (not saved) or forceReread:
    eig=np.array([readEigFile(systemLabel+str(i)+'td.EIG',sep=False) for i in steps])
    np.save('TDEIG',eig)
  
  return time, eig
  
#------------------------------------------------------------------- 
def getVASPBands():  
  f=open('EIGENVAL','r')
  for i in range(5):
  	f.readline()
  
  x,nk,nb=(int(i) for i in f.readline().split())
  #print nk,nb
  bands=np.zeros([nk,nb])
  kpoints=np.zeros([nk,3])
  x=range(nk)
  for i in range(nk):
  	f.readline()
  	coor=f.readline().split()
  	for k in range(3):
  		kpoints[i,k]=float(coor[k])	
  	for j in range(nb):
  		bands[i,j]=float(f.readline().split()[1])
  return x, bands

#------------------------------------------------------------------- 
def getBands(bandType='e'):
  """
  return the bands read from systemLabel.bands as the dimension numKpoint x numBand x numSpin
  returns: X, Ek, xticks, xticklabels
  """  
  import math
  options = tdapOptions()
  systemLabel = options.label
  band_file = open(systemLabel+'.bands')
  if bandType == 'e':
    FermiEnergy = float(band_file.readline().split()[0])
  else:
    band_file.readline()
    FermiEnergy = 0.0
  
  kmin,kmax=(float(i) for i in band_file.readline().split())
  emin,emax=(float(i) for i in band_file.readline().split())
  
  numBand,numSpin,numKpoint=(int(i) for i in band_file.readline().split())
  Ek = np.zeros([numKpoint,numBand,numSpin])
  X = np.zeros(numKpoint)
  
  for kpt in range(numKpoint):
    eigen = []
    for ispin in range(numSpin):
      for iband in range(int(math.ceil(numBand/10.0))):
        line = band_file.readline()
        eigen.extend([float(i) for i in line.split()])  
    X[kpt]=eigen[0]
    eigen.pop()
    Ek[kpt,:,:] = np.array(eigen).reshape(numBand,numSpin) - FermiEnergy #.reshape(1,-1) - EFermi
  
  numSpecK = int(band_file.readline().split()[0])
  xticks = []
  xticklabels =[]
  for i in range(numSpecK):
    line = band_file.readline()
    xticks.append(float(line.split()[0]))
    xticklabels.append(r'$'+line.split()[1][1:-1]+'$')
  
  return X, Ek, xticks, xticklabels
#------------------------------------------------------------------- 
def getDOS(filename = 'siesta.dos'):
  """
  return the dos read from filename as the dimension numEnergyPoint
  returns: energy, up, down, total
  """  
  dosfile = open(filename)
  for i in range(5):
    dosfile.readline()
    
  NumEinter=int(dosfile.readline().split()[2])  
  for i in range(3):
    dosfile.readline()  
  FermiEnergy=float(dosfile.readline().split()[3])
  
  for i in range(3):
    dosfile.readline()
    
  data = []
  for i in range(NumEinter):
    data.append([float(j) for j in dosfile.readline().split()])
   
  data = np.array(data) 
  energy = np.array(data[:,0]) - FermiEnergy
  
  up = np.array(data[:,1])
  down = np.array(data[:,2])
  total = np.array(data[:,3])
  
  return energy, up, down, total
  
#-------------------------------------------------------------------

  