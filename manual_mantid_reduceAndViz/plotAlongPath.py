import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from mantid.simpleapi import CreateMDHistoWorkspace
from mantid import plots

import sys
sys.path.append('/SNS/groups/dgs/DGS_SC_scripts')
from reduce_data_to_MDE import *
from slice_utils import *
import define_data
import define_slices
import os
import imp

# Manually plot the path for the FeSi data
a=4.489
#Gamma is (0,0,0)
#X is (0.5,0,0)
#M is (0.5,0.5,0)
#gamma_a=np.pi/d
#gamma_m=2.*np.pi/np.sqrt(3.)/a
#m_k=2.*np.pi/3/a
#gamma_k=4.*np.pi/3/a
m_gamma = np.sqrt(2)
gamma_r = np.sqrt(3)
r_x = np.sqrt(2)

plt.rcParams['figure.dpi'] = 300
plt.rcParams['savefig.dpi'] = 300
plt.rcParams.update({'font.size': 15})

fig=plt.figure(figsize=(12,5))
gs = GridSpec(1, 5,width_ratios=[1,1,m_gamma,gamma_r,r_x],wspace=0)

ax1 = plt.subplot(gs[0],projection='mantid')
ax2 = plt.subplot(gs[1],sharey=ax1,projection='mantid')
ax3 = plt.subplot(gs[2],sharey=ax1,projection='mantid')
ax4 = plt.subplot(gs[3],sharey=ax1,projection='mantid')
ax5 = plt.subplot(gs[4],sharey=ax1,projection='mantid')

colormesh_pars={}
colormesh_pars['norm']=LogNorm(vmin=1e-5,vmax=1e-2)

mtd.importAll()
ax1.pcolormesh(FeSi_300K_Ei70meV_GAMMA2X, **colormesh_pars)
ax2.pcolormesh(FeSi_300K_Ei70meV_X2M, **colormesh_pars)
ax3.pcolormesh(FeSi_300K_Ei70meV_M2GAMMA, **colormesh_pars)
ax4.pcolormesh(FeSi_300K_Ei70meV_GAMMA2R, **colormesh_pars)
ax5.pcolormesh(FeSi_300K_Ei70meV_R2X, **colormesh_pars)

#Adjust plotting parameters

ax1.set_ylabel('E (meV)')
ax1.set_xlabel('')
#ax1.set_xlim(0,1./3)
ax1.set_ylim(-10.,30.)
#ax1.set_title(r'$[\epsilon,\epsilon,0], 0 \leq \epsilon \leq 1/3$')
ax1.set_xticks([0,0.5])
ax1.set_xticklabels(['$\Gamma$','$X$'])
#ax1.spines['right'].set_visible(False)
ax1.tick_params(direction='in')

ax2.get_yaxis().set_visible(False)
#ax2.set_xlim(1./3,1./2)
ax2.set_xlabel('')
#ax2.set_title(r'$[1/3+\epsilon,1/3-2\epsilon,0], 0 \leq \epsilon \leq 1/6$')
ax2.set_xticks([0.5])
ax2.set_xticklabels(['$M$'])
#ax2.spines['left'].set_visible(False)
ax2.tick_params(direction='in')

ax3.get_yaxis().set_visible(False)
ax3.set_xlim(0.5,0)
ax3.set_xlabel('')
#ax2.set_title(r'$[1/3+\epsilon,1/3-2\epsilon,0], 0 \leq \epsilon \leq 1/6$')
ax3.set_xticks([0])
ax3.set_xticklabels(['$\Gamma$'])
#ax2.spines['left'].set_visible(False)
ax3.tick_params(direction='in')

ax4.get_yaxis().set_visible(False)
#ax2.set_xlim(1./3,1./2)
ax4.set_xlabel('')
#ax2.set_title(r'$[1/3+\epsilon,1/3-2\epsilon,0], 0 \leq \epsilon \leq 1/6$')
ax4.set_xticks([0.5])
ax4.set_xticklabels(['$R$'])
#ax2.spines['left'].set_visible(False)
ax4.tick_params(direction='in')

ax5.get_yaxis().set_visible(False)
#ax2.set_xlim(1./3,1./2)
ax5.set_xlabel('')
#ax2.set_title(r'$[1/3+\epsilon,1/3-2\epsilon,0], 0 \leq \epsilon \leq 1/6$')
ax5.set_xticks([0.5])
ax5.set_xticklabels(['$X$'])
#ax2.spines['left'].set_visible(False)
ax5.tick_params(direction='in')

#mid = (fig.subplotpars.right + fig.subplotpars.left)/2
#plt.suptitle('FeSi S(Q,E) along high symmetry path in 1st BZ',x=mid)

fig.show()


