#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 24 08:42:57 2018

@author: clian
"""

from ase.io import read

from py3ramids.io.output import writeQE

atoms = read('result',format='espresso-out')
writeQE('struct.part',atoms)
print(atoms)

