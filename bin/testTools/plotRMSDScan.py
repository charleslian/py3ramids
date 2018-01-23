import numpy as np
from matplotlib import pyplot as plt
import pyramids.io.result as dp
from pyramids.plot.PlotUtility import scanFolder
import pyramids.plot.setting as ma
import pyramids.process.struct as pps

#--------------------------------------------------------------------------------------------
def action(index,folder):
  dataCurFolder = []
  dataCurFolder.append([index, folder])
  dp.getTrajactory()
  dataCurFolder.append(pps.calculateRMSD())
  dataCurFolder.append(dp.getEnergyTemperaturePressure())
  return dataCurFolder

M = 28.08
R = 8.3144598
#--------------------------------------------------------------------------------------------
fig, axs = plt.subplots(2,1,sharex=True,sharey=False,figsize=(6,8))#
data = scanFolder(action,folders=['0.0025','0.0035','0.0045','0.0055','0.0065'])
c = ma.getColors(len(data))


melt_v = []
for line in data:
  index, folder = line[0]
  ax = axs[0]
  ax.plot(line[1][0], line[1][1], lw=3, label=folder, c=c[index])
  selectedstep = 5000
  t, v = float(folder)/2.0, line[1][1][selectedstep]/line[1][0][selectedstep]*1E3
  #melt_v.append(t, v)
  if index == 0:
      
      vt = np.sqrt(3*R*300/(M*1E-3))*1E-2
      ax.plot(line[1][0],line[1][0]*vt/1E3)
      print t, line[1][0][selectedstep], v, vt, v - vt
  kargs=ma.getPropertyFromPosition(xlabel='Time (fs)', ylabel=r'$\langle u \rangle^\frac{1}{2}$ ($\AA$)', 
                                   title='RMSD')
  ma.setProperty(ax,**kargs)
  #--------------------------------------------------------------------------------------------
  ax = axs[1]
  ax.plot(line[2][0], line[2][1], lw=3,c=c[index])
  kargs=ma.getPropertyFromPosition(xlabel='Time (fs)', ylabel=r'T (K)', 
                                   title='Temperature')
  ma.setProperty(ax,**kargs)
  
  
plt.tight_layout()
SaveName = __file__.split('/')[-1].split('.')[0]
if True:
  for save_type in ['.pdf','.png']:
    filename = SaveName + save_type
    #plt.savefig(filename,dpi=600)