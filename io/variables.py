"""
This module include the tdapOption class, 
which is initializied by reading the input.fdf file,
to get the useful setups, such as the system label, 
time step of electronic and ionic motions, spin status,
for the following data process. 
"""

__author__ = "Chao Lian <charleslian@126.com>"
__date__ = "Dec 5 2017"

__version__ = "3.0"
__credits__ = "Chao Lian initial and maintain the codes"

import numpy as np
from ase.io.trajectory import Trajectory
#import os

class Variables(object):
  def __init__(self):
    self.label  = ''
    
    self.numAtom = 1
    self.numElect = 1
    self.numKpt = 1
    self.numBnd = 1
    
    self.mdTimeStep = 0.0, 'fs'
    self.tdTimeStep = 0.0, 'fs'
    self.mdFinalStep = 1
    self.tdFinalStep = 1
    
    self.eigStartStep = 1
    self.eigSkipStep = self.mdFinalStep
    
    self.spinPol = False 
    self.curStep = 1      # run to this step
    self.eField = np.zeros([1,3])
    self.aField = np.zeros([1,3])
    self.kpts = np.zeros([1,4])
    self.traj = None      # ionic trajectory
    self.carrier = None   # excited electron and hole
    self.energy = None    # energies
    self.time = [0.0]
    
  def output(self):
    """
    Print the obtained values for check and test
    """
    
    print("System label is %s"%self.label)
    
    print("Total number of atoms is %i"%self.numAtom)
    print("Total number of electrons is %i"%self.numElect)
    print("Total number of k points is %i"%self.numKpt)
    print("Total number of bands is %i"%self.numBnd)
    
    print("Length of MD time step  %f %s"%(self.mdTimeStep[0],self.mdTimeStep[1]))
    print("Length of TD time step  %f %s"%(self.tdTimeStep[0],self.tdTimeStep[1]))
    print("Final MD time step is %i"%self.mdFinalStep)
    print("Final TD time step is %i"%self.tdFinalStep)
    
    print("System is spin polarized",self.spinPol)
    print("Current step %i, t = %f fs"%(self.curStep, self.time[-1]))
    
    #print("Kpoints are",self.kpts)
    

#-------------------------------------------------------------------
if __name__== '__main__':
  options = Variables()
  options.output()
