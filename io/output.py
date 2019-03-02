# -*- coding: utf-8 -*-
"""
Created on Wed Jul 27 00:10:57 2016

@author: clian
"""

listPAW = [
    "Ac.pbe-spfn-kjpaw_psl.1.0.0.UPF",
    "Ag.pbe-n-kjpaw_psl.1.0.0.UPF",
    "Ag.pbe-spn-kjpaw_psl.1.0.0.UPF",
    "Al.pbe-n-kjpaw_psl.1.0.0.UPF",
    "Al.pbe-nl-kjpaw_psl.1.0.0.UPF",
    "Am.pbe-spfn-kjpaw_psl.1.0.0.UPF",
    "Ar.pbe-n-kjpaw_psl.1.0.0.UPF",
    "Ar.pbe-nl-kjpaw_psl.1.0.0.UPF",
    "As.pbe-dn-kjpaw_psl.1.0.0.UPF",
    "As.pbe-n-kjpaw_psl.1.0.0.UPF",
    "At.pbe-dn-kjpaw_psl.1.0.0.UPF",
    "Au.pbe-n-kjpaw_psl.1.0.0.UPF",
    "Au.pbe-spfn-kjpaw_psl.1.0.0.UPF",
    "Au.pbe-spn-kjpaw_psl.1.0.0.UPF",
    "Ba.pbe-spn-kjpaw_psl.1.0.0.UPF",
    "Be.pbe-n-kjpaw_psl.1.0.0.UPF",
    "Be.pbe-s-kjpaw_psl.1.0.0.UPF",
    "Be.pbe-sl-kjpaw_psl.1.0.0.UPF",
    "Bi.pbe-dn-kjpaw_psl.1.0.0.UPF",
    "B.pbe-n-kjpaw_psl.1.0.0.UPF",
    "Br.pbe-dn-kjpaw_psl.1.0.0.UPF",
    "Br.pbe-n-kjpaw_psl.1.0.0.UPF",
    "Ca.pbe-spn-kjpaw_psl.1.0.0.UPF",
    "Cd.pbe-n-kjpaw_psl.1.0.0.UPF",
    "Cd.pbe-spn-kjpaw_psl.1.0.0.UPF",
    "Ce.pbe-spdfn-kjpaw_psl.1.0.0.UPF",
    "Ce.pbe-spdn-kjpaw_psl.1.0.0.UPF",
    "Cl.pbe-n-kjpaw_psl.1.0.0.UPF",
    "Cl.pbe-nl-kjpaw_psl.1.0.0.UPF",
    "Co.pbe-n-kjpaw_psl.1.0.0.UPF",
    "Co.pbe-spn-kjpaw_psl.1.0.0.UPF",
    "C.pbe-n-kjpaw_psl.1.0.0.UPF",
    "Cr.pbe-spn-kjpaw_psl.1.0.0.UPF",
    "Cs.pbe-spn-kjpaw_psl.1.0.0.UPF",
    "Cs.pbe-spnl-kjpaw_psl.1.0.0.UPF",
    "Cu.pbe-dn-kjpaw_psl.1.0.0.UPF",
    "Cu.pbe-spn-kjpaw_psl.1.0.0.UPF",
    "Dy.pbe-spdfn-kjpaw_psl.1.0.0.UPF",
    "Dy.pbe-spdn-kjpaw_psl.1.0.0.UPF",
    "Er.pbe-spdfn-kjpaw_psl.1.0.0.UPF",
    "Er.pbe-spdn-kjpaw_psl.1.0.0.UPF",
    "Eu.pbe-spdfn-kjpaw_psl.1.0.0.UPF",
    "Eu.pbe-spdn-kjpaw_psl.1.0.0.UPF",
    "Eu.pbe-spn-kjpaw_psl.1.0.0.UPF",
    "Fe.pbe-n-kjpaw_psl.1.0.0.UPF",
    "Fe.pbe-spn-kjpaw_psl.1.0.0.UPF",
    "F.pbe-n-kjpaw_psl.1.0.0.UPF",
    "Fr.pbe-spdn-kjpaw_psl.1.0.0.UPF",
    "Ga.pbe-dn-kjpaw_psl.1.0.0.UPF",
    "Ga.pbe-dnl-kjpaw_psl.1.0.0.UPF",
    "Gd.pbe-spdfn-kjpaw_psl.1.0.0.UPF",
    "Gd.pbe-spdn-kjpaw_psl.1.0.0.UPF",
    "Ge.pbe-dn-kjpaw_psl.1.0.0.UPF",
    "Ge.pbe-n-kjpaw_psl.1.0.0.UPF",
    "He.pbe-kjpaw_psl.1.0.0.UPF",
    "Hf.pbe-spdfn-kjpaw_psl.1.0.0.UPF",
    "Hf.pbe-spfn-kjpaw_psl.1.0.0.UPF",
    "Hf.pbe-spn-kjpaw_psl.1.0.0.UPF",
    "Hg.pbe-n-kjpaw_psl.1.0.0.UPF",
    "Hg.pbe-spn-kjpaw_psl.1.0.0.UPF",
    "Ho.pbe-spdfn-kjpaw_psl.1.0.0.UPF",
    "Ho.pbe-spdn-kjpaw_psl.1.0.0.UPF",
    "H.pbe-kjpaw_psl.1.0.0.UPF",
    "In.pbe-dn-kjpaw_psl.1.0.0.UPF",
    "I.pbe-dn-kjpaw_psl.1.0.0.UPF",
    "I.pbe-n-kjpaw_psl.1.0.0.UPF",
    "Ir.pbe-n-kjpaw_psl.1.0.0.UPF",
    "Ir.pbe-spfn-kjpaw_psl.1.0.0.UPF",
    "Ir.pbe-spn-kjpaw_psl.1.0.0.UPF",
    "K.pbe-spn-kjpaw_psl.1.0.0.UPF",
    "Kr.pbe-dn-kjpaw_psl.1.0.0.UPF",
    "La.pbe-spdfn-kjpaw_psl.1.0.0.UPF",
    "La.pbe-spfn-kjpaw_psl.1.0.0.UPF",
    "Li.pbe-s-kjpaw_psl.1.0.0.UPF",
    "Li.pbe-sl-kjpaw_psl.1.0.0.UPF",
    "Lu.pbe-spdfn-kjpaw_psl.1.0.0.UPF",
    "Lu.pbe-spdn-kjpaw_psl.1.0.0.UPF",
    "Mg.pbe-spn-kjpaw_psl.1.0.0.UPF",
    "Mg.pbe-spnl-kjpaw_psl.1.0.0.UPF",
    "Mn.pbe-spn-kjpaw_psl.1.0.0.UPF",
    "Mo.pbe-spn-kjpaw_psl.1.0.0.UPF",
    "Na.pbe-spn-kjpaw_psl.1.0.0.UPF",
    "Na.pbe-spnl-kjpaw_psl.1.0.0.UPF",
    "Nb.pbe-spn-kjpaw_psl.1.0.0.UPF",
    "Nd.pbe-spdfn-kjpaw_psl.1.0.0.UPF",
    "Nd.pbe-spdn-kjpaw_psl.1.0.0.UPF",
    "Ne.pbe-n-kjpaw_psl.1.0.0.UPF",
    "Ni.pbe-n-kjpaw_psl.1.0.0.UPF",
    "Ni.pbe-spn-kjpaw_psl.1.0.0.UPF",
    "N.pbe-n-kjpaw_psl.1.0.0.UPF",
    "Np.pbe-spfn-kjpaw_psl.1.0.0.UPF",
    "O.pbe-n-kjpaw_psl.1.0.0.UPF",
    "O.pbe-nl-kjpaw_psl.1.0.0.UPF",
    "Os.pbe-spfn-kjpaw_psl.1.0.0.UPF",
    "Os.pbe-spn-kjpaw_psl.1.0.0.UPF",
    "Pa.pbe-spfn-kjpaw_psl.1.0.0.UPF",
    "Pb.pbe-dn-kjpaw_psl.1.0.0.UPF",
    "Pd.pbe-n-kjpaw_psl.1.0.0.UPF",
    "Pd.pbe-spn-kjpaw_psl.1.0.0.UPF",
    "Pm.pbe-spdfn-kjpaw_psl.1.0.0.UPF",
    "Pm.pbe-spdn-kjpaw_psl.1.0.0.UPF",
    "Po.pbe-dn-kjpaw_psl.1.0.0.UPF",
    "P.pbe-n-kjpaw_psl.1.0.0.UPF",
    "P.pbe-nl-kjpaw_psl.1.0.0.UPF",
    "Pr.pbe-spdfn-kjpaw_psl.1.0.0.UPF",
    "Pr.pbe-spdn-kjpaw_psl.1.0.0.UPF",
    "Pt.pbe-n-kjpaw_psl.1.0.0.UPF",
    "Pt.pbe-spfn-kjpaw_psl.1.0.0.UPF",
    "Pt.pbe-spn-kjpaw_psl.1.0.0.UPF",
    "Pu.pbe-spfn-kjpaw_psl.1.0.0.UPF",
    "Ra.pbe-spdn-kjpaw_psl.1.0.0.UPF",
    "Rb.pbe-spn-kjpaw_psl.1.0.0.UPF",
    "Re.pbe-spfn-kjpaw_psl.1.0.0.UPF",
    "Re.pbe-spn-kjpaw_psl.1.0.0.UPF",
    "Rh.pbe-spn-kjpaw_psl.1.0.0.UPF",
    "Rn.pbe-dn-kjpaw_psl.1.0.0.UPF",
    "Ru.pbe-spn-kjpaw_psl.1.0.0.UPF",
    "Sb.pbe-dn-kjpaw_psl.1.0.0.UPF",
    "Sb.pbe-n-kjpaw_psl.1.0.0.UPF",
    "Sc.pbe-spn-kjpaw_psl.1.0.0.UPF",
    "Se.pbe-dn-kjpaw_psl.1.0.0.UPF",
    "Se.pbe-n-kjpaw_psl.1.0.0.UPF",
    "Si.pbe-n-kjpaw_psl.1.0.0.UPF",
    "Si.pbe-nl-kjpaw_psl.1.0.0.UPF",
    "Sm.pbe-spdfn-kjpaw_psl.1.0.0.UPF",
    "Sm.pbe-spdn-kjpaw_psl.1.0.0.UPF",
    "Sn.pbe-dn-kjpaw_psl.1.0.0.UPF",
    "S.pbe-n-kjpaw_psl.1.0.0.UPF",
    "S.pbe-nl-kjpaw_psl.1.0.0.UPF",
    "Sr.pbe-spn-kjpaw_psl.1.0.0.UPF",
    "Ta.pbe-spfn-kjpaw_psl.1.0.0.UPF",
    "Ta.pbe-spn-kjpaw_psl.1.0.0.UPF",
    "Tb.pbe-spdfn-kjpaw_psl.1.0.0.UPF",
    "Tb.pbe-spdn-kjpaw_psl.1.0.0.UPF",
    "Tc.pbe-spn-kjpaw_psl.1.0.0.UPF",
    "Te.pbe-dn-kjpaw_psl.1.0.0.UPF",
    "Te.pbe-n-kjpaw_psl.1.0.0.UPF",
    "Th.pbe-spfn-kjpaw_psl.1.0.0.UPF",
    "Ti.pbe-spn-kjpaw_psl.1.0.0.UPF",
    "Tl.pbe-dn-kjpaw_psl.1.0.0.UPF",
    "Tm.pbe-spdfn-kjpaw_psl.1.0.0.UPF",
    "Tm.pbe-spdn-kjpaw_psl.1.0.0.UPF",
    "U.pbe-spfn-kjpaw_psl.1.0.0.UPF",
    "V.pbe-spn-kjpaw_psl.1.0.0.UPF",
    "V.pbe-spnl-kjpaw_psl.1.0.0.UPF",
    "W.pbe-spfn-kjpaw_psl.1.0.0.UPF",
    "W.pbe-spn-kjpaw_psl.1.0.0.UPF",
    "Xe.pbe-dn-kjpaw_psl.1.0.0.UPF",
    "Yb.pbe-spdfn-kjpaw_psl.1.0.0.UPF",
    "Yb.pbe-spdn-kjpaw_psl.1.0.0.UPF",
    "Yb.pbe-spn-kjpaw_psl.1.0.0.UPF",
    "Y.pbe-spn-kjpaw_psl.1.0.0.UPF",
    "Zn.pbe-dn-kjpaw_psl.1.0.0.UPF",
    "Zn.pbe-spn-kjpaw_psl.1.0.0.UPF",
    "Zr.pbe-spn-kjpaw_psl.1.0.0.UPF",
]
def writeGPAWdipole(filename,efield,time,dipole):
  import numpy as np
  f = open(filename,'w')
  f.write('# Kick = '+str(efield)+'\n')
  f.write('#            time            norm                    dmx                    dmy                    dmz \n')
  for t ,dm in zip(time,dipole):
    norm = np.linalg.norm(dm)
    line = ('%20.8lf %20.8le %22.12le %22.12le %22.12le \n'
                    % (t, norm, dm[0], dm[1], dm[2]))
    f.write(line)
    
  # Kick = [    0.000000000000e+00,     1.000000000000e-05,     0.000000000000e+00]

def writeKLines(filename, points):
  lines = '%block BandLines \n'
  N = 50
  for index, kpt in enumerate(points):
    #print kpt
    out = (N*index+1, kpt[0][0], kpt[0][1], kpt[0][2], kpt[1])
    lines += '%i %3.2f %3.2f %3.2f %s \n'% out
  lines += '%endblock BandLines \n'  
  fdf = open(filename,'w')
  fdf.writelines(lines)
    
  
def writeSiesta(filename,atoms):
  NumberOfAtoms=atoms.get_number_of_atoms()
  elements=set(zip(atoms.get_chemical_symbols(), atoms.get_atomic_numbers()))
  #print list(elements).sort(key=) 
  NumberOfSpecies=len(elements)
  
  cell=atoms.get_cell()
  
  f=open(filename,'w')
  f.write("AtomicCoordinatesFormat  Ang\n")
  f.write("LatticeConstant  1.0  Ang\n\n")
  
  f.write("NumberOfAtoms    " + str(NumberOfAtoms)   + "\n")
  f.write("NumberOfSpecies  " + str(NumberOfSpecies) + "\n\n")
  
  f.write("%block LatticeVectors\n")
  lines= '' 
  for a in cell:
    lines += "  %21.16f %21.16f %21.16f\n" % tuple(a)
  f.write(lines)
  f.write("%endblock LatticeVectors\n\n")
  
  
  f.write("%block ChemicalSpeciesLabel\n")
  lines= ''
  
  element_index=dict()
  for i,element in enumerate(elements):
    element_index[element[0]]=i+1
    lines+= "%4d %5d %5s\n" % (i+1, element[1], element[0])
  
  #print element_index
  f.write(lines)
  f.write("%endblock ChemicalSpeciesLabel\n\n")
  
  f.write("%block AtomicCoordinatesAndAtomicSpecies\n")
  lines= ''
  for a in zip(atoms.get_positions(),atoms.get_chemical_symbols()):
    lines += " %21.16f %21.16f %21.16f " % tuple(a[0])
    lines += "%2d\n" % element_index[a[1]]
    
  f.write(lines)
  f.write("%endblock AtomicCoordinatesAndAtomicSpecies\n\n")
  
def splitMDCAR():
  """
  split the siesta.MD_CAR file to POSCAR file per step
  """
  import os
  systemLabel = 'siesta'
  NumBlocks=int(os.popen('grep -i '+systemLabel+' '+systemLabel+'.MD_CAR | wc -l').readline().split()[0])
  position_file = open('siesta'+'.MD_CAR')
  atomNumList = [int(i) for i in os.popen('head -6 siesta.MD_CAR |tail -1').readline().split()]
  numAtomPositionLine = sum(atomNumList)
  totalNumLine = numAtomPositionLine + 7
  context = position_file.readlines()
  for index in range(NumBlocks):
    output=context[index*totalNumLine:(index+1)*totalNumLine]
    poscarFileName = "POSCAR"+str(index)
    poscarFile=open(poscarFileName,'w')
    poscarFile.writelines(output)
    
def writeQE(filename,atoms):
  NumberOfAtoms=atoms.get_number_of_atoms()
  elements=set(zip(atoms.get_chemical_symbols(), 
                   atoms.get_atomic_numbers(), atoms.get_masses()))
  #print list(elements).sort(key=) 
  NumberOfSpecies=len(elements)
  cell=atoms.get_cell()
  #print atoms.get_masses()
  
  f=open(filename,'w')

  f.write('    ibrav = 0, ')
  f.write('  nat = ' + str(NumberOfAtoms) + ",")
  f.write('  ntyp = ' + str(NumberOfSpecies) + ",\n")
  
  f.write("CELL_PARAMETERS angstrom\n")
  lines= '' 
  for a in cell:
    lines += "  %21.16f %21.16f %21.16f\n" % tuple(a)
  f.write(lines)
  
  #f.write("%block ChemicalSpeciesLabel\n")

  #print element_index
  
  #f.write("%endblock ChemicalSpeciesLabel\n\n")
  
  
  lines= 'ATOMIC_SPECIES\n'
  for i,element in enumerate(elements):
    psedo = [ele for ele in listPAW if ele.split('.')[0]==element[0]][0]
    lines+= "%5s %21.16f %s\n" % (element[0], element[2], psedo)
  f.write(lines)
  
  
  f.write("ATOMIC_POSITIONS angstrom\n")
  lines= ''
  for a in zip(atoms.get_positions(),atoms.get_chemical_symbols()):
    lines +=  a[1]
    lines += " %21.16f %21.16f %21.16f \n" % tuple(a[0])
    
  f.write(lines)
  

  
  f.write("\n\n")  
  
  
