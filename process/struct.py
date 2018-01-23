import numpy as np
def calculateRMSD(traj, selectedStep=None, atomsOrigin=None, init=0, selectedAtoms=None):
  """ 
  return the radius mean square displacements of the selected steps 
  compared with the init step
  """
  #from ase.io.trajectory import Trajectory
  if atomsOrigin is None:
    atomsOrigin = traj[init]
    #print 'no'
  if selectedStep == None:
    selectedStep = range(len(traj))

  distance = np.array([np.mean(calculateDisplacement(traj[step],atomsOrigin,selectedAtoms)**2)**0.5
              for step in selectedStep])
  return distance

#-------------------------------------------------------------------
def calculateDisplacement(atoms,atomsOrigin,selectedAtoms=None):
  """ 
  return the displacements as a dimension of Natoms
  """
  if selectedAtoms == None:
    selectedAtoms = range(atomsOrigin.get_number_of_atoms())
    
  return np.array([
          np.min([np.linalg.norm(atoms.get_positions()[index] - 
                  atomsOrigin.get_positions()[index] + np.dot([i,j,k],atoms.get_cell())) 
                  for i in range(-1,2) 
                  for j in range(-1,2) 
                  for k in range(-1,2)
                  ])
          for index in selectedAtoms])