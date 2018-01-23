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

class TdpwVarible(Variables):
  def __init__(self,inputFile='input.in', outputFile = 'result'):
    self.inputFile  = inputFile
    self.outputFile = outputFile
    self.eigStartStep = 2
    #self.options    = self.readFdf(self.inputFile)
    self.label      = 'pwscf'
    self.mdTimeStep = self._readValue('dt',float)*0.048378, 'fs'
    self.tdTimeStep = self._readValue('edt',float)*0.048378, 'fs'
    self.mdFinalStep = self._readValue('nstep',int)
    self.tdFinalStep = self._readValue('nstep',int)
    self.spinPol = False #self.__getLogical('spinpolarized',False)
    self.eigLengthStep = 1
    #self.laserParam = self.__getArray('tdlightenvelope',np.zeros(5))
    
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
    self.carrier = self.getCarrier()
    self.cartime = self.time
  
  def _readValue(self, tag, typeFunc):
    value = self._findLineContain(tag, self.inputFile)[1].split()[-1]
    if value[-1] == ',' : value = value[:-1]
    return typeFunc(value)
  
  def getEfield(self):
    with open('TDEFIELD') as file: 
      return np.array([[float(i) for i in line.split()] for line in file])*1E-5 * Rydberg/Bohr

  def getAfield(self):
    with open('TDAFIELD') as file: 
      return np.array([[float(i) for i in line.split()] for line in file])*1E-5 #* Rydberg/Bohr
      
  def getEnergy(self): 
    """
    FIXME: not read KS energy
    return the Temperature, KS Energy and Total Energy as the dimension of Nstep
    read from systemLabel.MDE
    returns the Temperature, KS Energy, Total Energy, Volume, Pressure
    """
    lines = self._findAllLineContain("!    total energy", self.outputFile)
    ksEnergy = np.array([float(line.split()[-2])*Rydberg for line in lines])
    
    lines = self._findAllLineContain('Ekin + Etot', self.outputFile)
    if lines[0].split()[-2] != 'NaN':
      totalEnergy = np.array([float(line.split()[-2])*Rydberg for line in lines])
    else:
      totalEnergy = ksEnergy
    
    
    return ksEnergy, totalEnergy
    
    
  def getTemperature(self): 
    """
    return the Temperature, KS Energy and Total Energy as the dimension of Nstep
    read from systemLabel.MDE
    returns the Temperature, KS Energy, Total Energy, Volume, Pressure
    """
    lines = self._findAllLineContain("temperature           =", self.outputFile)
    temperature = np.array([float(line.split()[-2]) for line in lines])
    return temperature
  
  def _findLineContain(self, phrase, file):
    with open(file) as f:
      for index, line in enumerate(f):
        if phrase in line: return index, line

  def _findAllLineContain(self, phrase, file):
    with open(file) as f:
      return [line for line in f if phrase in line]
        
  def getNumAtom(self):
    keyword = 'number of atoms/cell'
    index, line = self._findLineContain(keyword,self.outputFile)
    return int(line.split()[-1])

  def getNumKpt(self):
    keyword = 'number of k points='
    index, line = self._findLineContain(keyword,self.outputFile)
    return int(line.split()[-1])
  
  def getNumElect(self):
    keyword = 'number of electrons'
    index, line = self._findLineContain(keyword,self.outputFile)
    return float(line.split()[-1])
  
  def getNumBnd(self):
    keyword = 'number of Kohn-Sham states'
    index, line = self._findLineContain(keyword,self.outputFile)
    return int(line.split()[-1])
  
  def getKpts(self):
    keyword = 'number of k points='
    index, line = self._findLineContain(keyword,self.outputFile)
    start = index+1
    with open(self.outputFile) as f:
      context = [temp.split() for i, temp in enumerate(f) if start+self.numKpt>=i>start]
        
        
   #print(context)
    kpts = np.array([[float(i) for i in (line[4], line[5], line[6][:-2], line[9])] for line in context])
      
    return kpts

  def getTrajactory(self, fileobj='result'):
    from ase.atoms import Atoms, Atom
    from ase import units
    from ase.calculators.singlepoint import SinglePointCalculator
    
    """Reads quantum espresso output text files."""
    
    fileobj = open(fileobj, 'rU')
    lines = fileobj.readlines()
    #print lines
    images = []
    
    for number, line in enumerate(lines):
      if 'number of atoms/cell' in line:
        numAtom = int(line.split()[-1])
      if 'number of atomic types' in line:
        numElement = int(line.split()[-1])
      
      if "atomic species   valence    mass" in line:
        elementLines = lines[number+1:number+1+numElement]
        elementsInfo = [[(i) for i in l.split()] for l in elementLines]
  
        #print numAtom
      if 'lattice parameter (alat)' in line:
        line = line.strip().split('=')[1].strip().split()[0]
        lattice_parameter = float(line) * units.Bohr
      
      if "crystal axes: (cart. coord. in units of alat)" in line:
        ca_line_no = number
        cell = np.zeros((3, 3))
        for number, line in enumerate(lines[ca_line_no + 1: ca_line_no + 4]):
            line = line.split('=')[1].strip()[1:-1]
            values = [float(value) for value in line.split()]
            cell[number, 0] = values[0]
            cell[number, 1] = values[1]
            cell[number, 2] = values[2]
        cell *= lattice_parameter
      if "site n.     atom                  positions (alat units)" in line:
        initPositionLines = lines[number+1:number+1+numAtom]
        elements = np.array([l.split()[1] for l in initPositionLines])
        scaledPositions = np.array([[float(i) for i in l.split()[6:9]] 
                              for l in initPositionLines])
        atom = Atoms(symbols=elements, scaled_positions=scaledPositions, cell=cell, pbc=True)
        images.append(atom)                       
             
                     
      cellkey = 'CELL_PARAMETERS (angstrom)'
      posikey = 'ATOMIC_POSITIONS' 
      if posikey in line:
        positionLines = lines[number+1:number+1+numAtom]
        elements = np.array([l.split()[0]
                          for l in positionLines])
        positions = np.array([[float(i) for i in l.split()[1:]] 
                          for l in positionLines])
        if 'angstrom' in line:
          atom = Atoms(symbols=elements, positions=positions, cell=cell, pbc=True)
        elif 'crystal' in line:
          atom = Atoms(symbols=elements, scaled_positions=positions, cell=cell, pbc=True)
        images.append(atom)
      
      if "Forces acting on atoms (cartesian axes, Ry/au):" in line:
        #print number
        atom = images[-1]
        forceLines = lines[number+2:number+2+numAtom]
        forces = np.array([[float(i) for i in l.split()[-3:]] 
                          for l in forceLines])
        forces *= units.Ry / units.Bohr
        calc = SinglePointCalculator(atom, forces=forces)
        atom.set_calculator(calc)
        
    filename = 'Trajectory'
    from ase.io import write  
    write(filename,images,'traj')
    
    from ase.io.trajectory import Trajectory
    
    traj = Trajectory(filename)
    return traj
    
  
  def getCurrent(self):
    keyword = 'current is'
    lines = self._findAllLineContain(keyword, self.outputFile) + ['current is 0.0 0.0 0.0']
    current = np.array([[float(num) for num in line.split()[2:6]] for line in lines])
    #print(current)
    current[-1,:] = current[-2,:]  
    return current - current[0,:] 
  
  def getCarrier(self):
    carrier = self.readTDData(mode='norm')
    carrier -= carrier[0,:,:]
    nocc = int(self.numElect)//2
    #print(self.numElect)
    hotHole = carrier[:,:,:nocc].sum(axis=(1,2))
    hotElectron = carrier[:,:,nocc:].sum(axis=(1,2))
    return hotElectron, hotHole
    
    
  
  def readTDData(self, mode='value'):
    if mode == 'norm':
      filename = 'pwscf.norm.dat'
      weighted = True
    elif mode == 'value':
      filename = 'pwscf.value.dat'
      weighted = False
      
    f = open(filename)
    text = f.readlines()
    nbnd, nkstot = [int(i) for i in text[0].split()]
    kweight = [float(i) for i in text[1].split()]
    
    nstep = (len(text) - 2)//(nkstot+1)
    del text[1]
    del text[::(nkstot+1)]
    data = np.array([[float(i) for i in line.split()] for line in text])
    data = data[:nstep*nkstot].reshape([nstep, nkstot, nbnd])
    
    if weighted:
      for ib in range(data.shape[2]):
        data[:,:,ib] *= kweight
      
    return data

#-------------------------------------------------------------------
if __name__== '__main__':
  options = TdpwVarible()
  options.output()
