#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 11 08:41:57 2019

@author: clian
"""

from ase.io import read
from py3ramids.io.output import writeQE

def generate(filename='input.in', filetype='espresso-in',superindex=(2,2,2)):
    atoms = read(filename,format=filetype)
    writeQE('supercell.in',atoms*superindex)


if __name__ == '__main__':
  import argparse
  
  class options(argparse.ArgumentParser):
    def __init__(self):
      super(options, self).__init__()
      self.add_argument('-s', '--superindex', type=int, nargs='+',  
                        help='cell index')
      self.args = self.parse_args() 
      
  options = options()
  args = options.args
  generate(superindex=tuple(args.superindex))
  #print()
#      self.add_argument('-f', '--frequency', nargs=None,  
#                        help='frequency as input, unit in 1/fs')
#      self.add_argument('-l', '--wavelength', nargs=None,  
#                        help='wavelength as input, unit in nm')
#      self.add_argument('-t', '--period', nargs=None,  
#                        help='period as input, unit in fs')                   
#      self.add_argument('-e', '--energy', nargs=None,  
#                        help='energy as input, unit in eV')                        
#      self.args = self.parse_args() 