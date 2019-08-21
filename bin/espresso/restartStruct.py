#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 24 08:42:57 2018

@author: clian
"""

from ase.io import read, write
from py3ramids.io.output import writeQE, readQE
from ase.visualize import view
import argparse
  
class options(argparse.ArgumentParser):
    def __init__(self):
        super(options, self).__init__()
        self.add_argument('-v', '--view', type=bool, nargs='?', default=False,
                          help='view the structure')
        self.add_argument('-p', '--print', type=bool, nargs='?', default=False,
                          help='print the structure')
        self.add_argument('-cart', '--cartesian', type=bool, nargs='?', default=True,
                          help='output as cartesian')
        self.add_argument('-i', '--input', type=str, nargs='?', default='result',
                          help='input file, default: "result"')  
        self.add_argument('-o', '--output', type=str, nargs='?', default='struct.part',
                          help='output file, default: "struct.part"')
        self.add_argument('-pp', '--pseudopotential', type=str, nargs='?', default='PAW',
                          help='pseudopotential type manual, paw, nc, default: "paw"')
        self.args = self.parse_args()
    
options = options()
args = options.args
#print(args)
#var = TdpwVarible()
if args.pseudopotential == 'manual':
    lines = open('input.in').readlines()
    for il, line in enumerate(lines):
        if 'ATOMIC_SPECIES' in line:
            start = il
        if 'ATOMIC_POSITIONS' in line:
            end = il
    args.pseudopotential = ''.join(lines[start+1:end])
    #print(args.pseudopotential)
            
            
#atoms = var.getTrajactory()[-1]
#print(readQE(args.input).__n)
atoms = list(readQE(args.input))[-1]
#atomList = [atoms for atoms in readQE(args.input)]
    #atoms = i
#atoms = read(args.input,format='espresso-out')
if args.view:
    view(atoms)
if args.print:
    print(atoms)
writeQE(args.output,atoms,cart=args.cartesian,pp=args.pseudopotential)
write('struct.xsf',atoms)
#print(atoms)

