"""
This module include the tdapOption class, 
which is initializied by reading the input.fdf file,
to get the useful setups, such as the system label, 
time step of electronic and ionic motions, spin status,
for the following data process. 
"""

__author__ = "Chao Lian <charleslian@126.com>"
__date__ = "Fri July 22 18:16:20 2016"

__version__ = "2.0"
__credits__ = "Chao Lian initial and maintain the codes"

import numpy as np

from py3ramids.io.variables import Variables
from ase.units import Rydberg, Bohr

class TdapVarible(Variables):
  def __init__(self,inputFile='input.fdf', outputFile = 'result'):
    self.inputFile  = inputFile
    self.outputFile = outputFile
    self.eigStartStep = 3
    self.options    = self.readFdf(self.inputFile)
    self.label      = self.__getString('systemlabel','siesta')
    self.mdTimeStep = self.__getFloatWithUnit('mdlengthtimestep',(0.0,'fs'))
    self.tdTimeStep = self.__getFloatWithUnit('tdlengthtimestep',(0.025,'fs'))
    self.mdFinalStep = self.__getInteger('mdfinaltimestep',1)
    self.tdFinalStep = self.__getInteger('tdfinaltimestep',self.mdFinalStep)
    self.spinPol = self.__getLogical('spinpolarized',False)
    self.eigLengthStep = self.__getInteger('tdwriteeigpairstep',self.mdFinalStep)
    self.laserParam = self.__getArray('tdlightenvelope',np.zeros(5))
    
    # Basic 
    self.numAtom = self.getNumAtom()
    self.numElect = self.getNumElect()
    self.numKpt = self.getNumKpt()
    self.numBnd = self.getNumBnd()
    
    # Arrays 
    self.eField = self.getEfield()
    self.aField = self.getAfield()
    self.energy = self.getEnergy()
    
    self.curStep = self.energy[0].shape[0] # Derivatives
    self.time = np.arange(self.curStep)*self.tdTimeStep[0] # Derivatives
    
    self.temperature = self.getTemperature()
    self.kpts = self.getKpts()
    self.traj = self.getTrajactory()
    self.current = self.getCurrent()
    self.cartime, self.carrier = self.getCarrier()
  
  def getEfield(self):
    with open('TDEFIELD') as file: 
      return np.array([[float(i) for i in line.split()] for line in file])*1E-5 * Rydberg/Bohr

  def getAfield(self):
    with open('TDAFIELD') as file: 
      return np.array([[float(i) for i in line.split()] for line in file])*1E-5 #* Rydberg/Bohr
      
  def getEnergy(self): 
    """
    return the Temperature, KS Energy and Total Energy as the dimension of Nstep
    read from systemLabel.MDE
    returns the Temperature, KS Energy, Total Energy, Volume, Pressure
    """
    start = 0
    with open(self.label+'.MDE') as file:      
      data = np.array([[float(num.replace('*','0')) for num in line.split()]
                     for line in file.readlines()[1:]])[start:,:]
    totEenergy = data[:,2]
    KSEnergy = data[:,3]
    return totEenergy, KSEnergy
  
  def getTemperature(self): 
    """
    return the Temperature, KS Energy and Total Energy as the dimension of Nstep
    read from systemLabel.MDE
    returns the Temperature, KS Energy, Total Energy, Volume, Pressure
    """
    start = 0
    with open(self.label+'.MDE') as file:      
      data = np.array([[float(num.replace('*','0')) for num in line.split()]
                     for line in file.readlines()[1:]])[start:,:]
    temperature = data[:,1]
    return temperature
  
  def _findLineContain(self, phrase, file):
    with open(file) as f:
      for line in f:
        if phrase in line: return line

  def _findAllLineContain(self, phrase, file):
    with open(file) as f:
      return [line for line in f if phrase in line]
        
  def getNumAtom(self):
    keyword = 'Number of atoms, orbitals, and projectors'
    line = self._findLineContain(keyword,self.outputFile)
    return int(line.split()[-3])

  def getNumKpt(self):
    keyword = 'Number of k-points'
    line = self._findLineContain(keyword,self.outputFile)
    return int(line.split()[-1])
  
  def getNumElect(self):
    keyword = 'Total number of electrons:'
    line = self._findLineContain(keyword,self.outputFile)
    return float(line.split()[-1])
  
  def getNumBnd(self):
    keyword = 'Number of atoms, orbitals, and projectors'
    line = self._findLineContain(keyword,self.outputFile)
    return int(line.split()[-2])
  
  def getKpts(self):
    with open(self.label+'.KP') as f:
      kpts = np.array([[float(value) for value in line.split()[1:]] 
                        for line in f.readlines()[1:]])
      kpts[:,:3] *= 1.0/Bohr
    return kpts
  
  def getTrajactory(self):
    atoms = self.xv_to_atoms(self.label+'.XV')
    atoms.pbc = [True,True,True]
    filename = 'Trajectory'
    atomsList = []
    with open(self.label+'.MD_CAR') as f:
      blockLine = self.numAtom + 7
      context = f.readlines()
      
      for step in range(self.curStep-1):
        output = context[step*blockLine:(step+1)*blockLine]
        coodinates=np.array([[float(value.replace('\n','')) 
                                for value in line.split()] 
                                for line in output[7:]])
        atomsCurrent = atoms.copy()
        atomsCurrent.set_scaled_positions(coodinates)
        atomsList.append(atomsCurrent)   
        
    from ase.io import write  
    write(filename,atomsList,'traj')
    
    from ase.io.trajectory import Trajectory
    
    traj = Trajectory(filename)
    return traj
  
  def xv_to_atoms(self, filename):
    """Create atoms object from xv file.

    Parameters:
        -filename : str. The filename of the '.XV' file.

    return : An Atoms object
    """
    from ase.atoms import Atoms
    with open(filename, 'r') as f:
        # Read cell vectors (lines 1-3)
        vectors = []
        for i in range(3):
            data = (f.readline()).split()
            vectors.append([float(data[j]) * Bohr for j in range(3)])

        # Read number of atoms (line 4)
        int(f.readline().split()[0])

        # Read remaining lines
        speciesnumber, atomnumbers, xyz, V = [], [], [], []
        for line in f.readlines():
            if len(line) > 5:  # Ignore blank lines
                data = line.split()
                speciesnumber.append(int(data[0]))
                atomnumbers.append(int(data[1]))
                xyz.append([float(data[2 + j]) * Bohr for j in range(3)])
                V.append([float(data[5 + j]) * Bohr for j in range(3)])
                
    vectors = np.array(vectors)
    atomnumbers = np.array(atomnumbers)
    xyz = np.array(xyz)
    atoms = Atoms(numbers=atomnumbers, positions=xyz, cell=vectors)
    return atoms
  
  def getCurrent(self):
    keyword = 'TDAP: Afield: Current ='
    linesApp = 2*['TDAP: Afield: Current =       -0.000       -0.000       -0.000 a.u.\n']
    lines = self._findAllLineContain(keyword, self.outputFile) + linesApp
    if len(lines) == 2:
       return np.zeros([self.time.shape[0], 3])
    #print(lines)
    current = np.array([[float(num) for num in line.split()[4:7]] for line in lines])
    #print(current)
    current[-1,:] = current[-2,:] = current[-3,:]
    return current - current[0,:]
  
  def readTDData(self, mode='value'):
    if mode == 'norm':
      appendix = 'q.EIG'
      sep = False
    elif mode == 'value':
      appendix = '.EIG'
      sep = False
    #nocc = int(self.numElect)//2
    #kweight = self.kpts[:,3]
    selectStep = self.getEIGSteps()
    print(selectStep)
    return np.array([self.readEigFile(self.label+str(step)+appendix,sep=sep) for step in selectStep])
    
  #-------------------------------------------------------------------    
  def readEigFile(self, filename = 'siesta.EIG', sep = False):
    """
    return the Eigenvalues read from systemLabel.EIG
    """
    import math
    eigFile = open(filename)
    line = eigFile.readline()
    EFermi = float(line.split()[0])
    line = eigFile.readline()
    nband, nspin, nkpt= [int(i) for i in line.split()[:3]]
    eigen = []
    for kpt in range(nkpt):
      for ispin in range(nspin):
        eigenPerKpt = []
        for iband in range(int(math.ceil(nband/10.0))):
          line = eigFile.readline()
          eigenPerKpt.extend([float(i) for i in line.split()])     
      eigenPerKpt.pop(0)
      eigen.append(eigenPerKpt)
    eigenvalues = np.array(eigen)
    #print(eigenvalues)
    if sep:
      return eigenvalues, EFermi
    else:
      return eigenvalues - EFermi
    
  def getEIGSteps(self):
    steps = []
    import os
    for i in os.listdir('.'):
        if i[:6] == 'siesta' and i[-5:] == 'q.EIG':
            steps.append(int(i[6:-5]))
  
    steplist = np.array(np.sort(steps),dtype=int)  
    return steplist
  
  def getCarrier(self, selectK=None, comp = False, ave=False):
    """
    """
    nocc = int(self.numElect)//2
    kweight = self.kpts[:,3]
    selectStep = self.getEIGSteps()
    
    if len(selectStep) == 0:
      return [], [[],[]]
    
    #print(selectStep)
    timestep = self.tdTimeStep[0]
    selectTime = selectStep*timestep
    
    
    hotE = np.zeros([len(selectStep)])
    hotH = np.zeros([len(selectStep)])
    for index,step in enumerate(selectStep):
      partition = self.readEigFile(self.label+str(step)+'q.EIG')
      for i in range(partition.shape[0]):
        partition[i,:] *= kweight[i]
      if selectK is not None: 
        hotE[index] = np.sum(partition[selectK,nocc:])
      else:
        hotH[index] = nocc*2 - np.sum(partition[:,:nocc])
        hotE[index] = np.sum(partition[:,nocc:])
      
    #numAtom = getNumOfAtoms()
    #hotH /= self.numAtom
    #hotE /= self.numAtom
    return selectTime, (hotE, hotH)
    
  def readFdf(self,inputFile='input.fdf'):
    """
    Read the fdf file input.fdf.
    The tags, including both one-line tags and block tags, 
    are parsed by a dictonary, of which the keys are the tags 
    and the values are a string of the input value.  
    """
    fdfFile = open(inputFile)
    inBlock = False
    options = {}
    
    for line in fdfFile.readlines():
      # skip blank lines, containing only the '\n'      
      if len(line) == 1:
        continue
      if line[0] == '#':
        continue
      line = line.lower().replace('_','')
      line = line.split('#')[0]
      lineList = line.split()
      
      # Read the Blocks
      if lineList[0] == '%block':
        inBlock = True
        blockName = lineList[1].replace('.','')
        blockValues = []
        continue
      elif lineList[0] == '%endblock':
        inBlock = False
        options[blockName] = blockValues
        continue
      if inBlock:
        blockValues.append(lineList)
        
      # Read the line mode parameters 
      else:
        tag = lineList[0].replace('.','')
        if len(lineList) == 1:
          options[tag] = 'True'
        elif len(lineList) == 2:
          options[tag] = lineList[1]
        elif len(lineList) == 3:
          options[tag] = (lineList[1],lineList[2])
    return options
  #def output(self): pass
      #output(super)
 
  def __getLogical(self,label,default):
    if label in self.options.keys():
      value = self.options[label]
      if 't' in value.lower():
        return True
      else:
        return False
    else:
      return default
 
  def __getString(self,label,default):
    if label in self.options.keys():
      return self.options[label]
    else:
      return default
      
  def __getInteger(self,label,default):
    if label in self.options.keys():
      return int(self.options[label])
    else:
      return default

  def __getFloat(self,label,default):
    if label in self.options.keys():
      return float(self.options[label].replace('f','e'))
    else:
      return default
      
  def __getFloatWithUnit(self,label,default):
    if label in self.options.keys():
      return float(self.options[label][0].replace('f','e')), self.options[label][1]
    else:
      return default
      
  def __getArray(self,label,default):
    if label in self.options.keys():
      value = np.array([[float(i) for i in line] for line in self.options[label]])
      return value
    else:
      return default

#-------------------------------------------------------------------
if __name__== '__main__':
  options = TdapVarible()
  options.output()
