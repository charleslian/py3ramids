#/usr/bin/python3

import py3ramids.io.output as tdio
from ase.io import write, read

import argparse
class options(argparse.ArgumentParser):
    def __init__(self):
      super(options, self).__init__()
      
      self.add_argument('input', metavar='x', type=str, nargs='?',  
                        help='input value.')
      self.add_argument('-f', '--frequency', nargs=None,  
                        help='frequency as input, unit in 1/fs')
      self.add_argument('-l', '--wavelength', nargs=None,  
                        help='wavelength as input, unit in nm')
      self.add_argument('-t', '--period', nargs=None,  
                        help='period as input, unit in fs')                   
      self.add_argument('-e', '--energy', nargs=None,  
                        help='energy as input, unit in eV')                        
      self.args = self.parse_args()
    
#options = options()
#args = options.args
#atoms = read(args.input,format='xyz')#xv_to_atoms('siesta.XV')
atoms = read('input.xyz',format='xyz')
print(atoms.get_positions())
import numpy as np
#atoms.rotate([1,0,0],-0.5*np.pi,rotate_cell=True)
#atoms.pbc = [True, True, True]
#ase.build.niggli_reduce(atoms)
tdio.writeSiesta('structure.fdf',atoms)
tdio.writeQE('structure.in',atoms)
#write('struture.pdb',atoms,format='pdb')
#write('POSCAR',atoms)

#view(atoms*[2,2,2])