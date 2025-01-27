# code repurposed from pathSQE plot_along_path_foldedBZ()
import matplotlib.pyplot as plt
from matplotlib.colors import SymLogNorm
from matplotlib.gridspec import GridSpec
from mantid import plots
from mantid.simpleapi import *



def IqE_to_SED(SED_ws,T,Ebins):
    SED_ws_shape = SED_ws.getSignalArray().shape
    k = 8.617e-2
    # Split the string to extract min, step, and max values
    min_val, step_val, max_val = map(float, Ebins.split(','))

    # Generate the vector using numpy's arange function
    E_vals = np.arange(min_val + step_val/2, max_val + step_val/2, step_val)
    E_terms = np.exp(E_vals/(k*T))
    factors = E_vals / (1 + (1/(E_terms-1)))
    reshaped_factors = np.tile(factors.reshape(1, 1, 1, len(factors)), (SED_ws_shape[0], SED_ws_shape[1], SED_ws_shape[2], 1))
    SED_ws.setSignalArray(SED_ws.getSignalArray() * reshaped_factors)
        
    return SED_ws



outputs_folder = '/SNS/CNCS/IPTS-31965/shared/Aiden/comprehensive/BZfold_pathSQE/300K_20meV_initial/'

prim2mantid_transMatrix = np.linalg.inv(np.array([[2/3, -1/3, -1/3], [1/3, 1/3, -2/3], [1/3, 1/3, 1/3]])) # bilbao space group R-3m (166) to hexagon>

user_defined_path = {'path':[ ('Gamma','T'), ('T','X'), ('X','L'), ('L','Gamma'), ('Gamma','X') ], 
                         # defined in terms of the primitive rhombohedral basis (i.e. not the basis being used in mantid which is often the conventional basis)
                         'point_coords':{'Gamma': [0.0, 0.0, 0.0], 'T': [0.5, 0.5, 0.5], 'X': [0 , 0.5, 0.5], 'L': [0, 0.5, 0] }}

T = 300
Ebins = '-15,0.2,15'



foldedSegNames = []
foldedSEDNames = []
ratios = []
for seg in user_defined_path['path']:
    name = '{}2{}'.format(seg[0],seg[1])
    foldedSegNames.append(name)
    LoadMD(Filename=outputs_folder+'{}_folded.nxs'.format(name), OutputWorkspace='{}'.format(name), LoadHistory=False)
    
    SED_ws = mtd['{}'.format(name)].clone()
    SED_ws = IqE_to_SED(SED_ws=SED_ws,T=T,Ebins=Ebins)
    SaveMD(SED_ws, Filename=outputs_folder+'{}_folded_SEDnew.nxs'.format(name))
    LoadMD(Filename=outputs_folder+'{}_folded_SEDnew.nxs'.format(name), OutputWorkspace='{}_SED'.format(name), LoadHistory=False)
    foldedSEDNames.append('{}_SED'.format(name))
    
    pt1_prim = np.array(user_defined_path['point_coords'][seg[0]])
    pt2_prim = np.array(user_defined_path['point_coords'][seg[1]])
    norm_prim = np.linalg.norm(pt2_prim-pt1_prim)

    pt1_conv = prim2mantid_transMatrix @ pt1_prim
    pt2_conv = prim2mantid_transMatrix @ pt2_prim
    norm_conv = np.linalg.norm(pt2_conv-pt1_conv)

    ratios.append(norm_prim)

print(foldedSegNames)

print(foldedSEDNames)

print('scaling ratios ', ratios)

lens = [0, 0.126454, 0.260268, 0.394082, 0.527896, 0.680331]
ratios = [lens[i+1]-lens[i] for i in range(len(lens)-1)]


plt.rcParams['figure.dpi'] = 300
plt.rcParams['savefig.dpi'] = 300
plt.rcParams.update({'font.size': 28})

colormesh_pars={}
colormesh_pars['norm']=SymLogNorm(linthresh=5e-3,vmin=5e-3,vmax=2e-1) # data (linthresh=1e-3,vmin=1e-3,vmax=2e-1)
colormesh_pars['cmap']='viridis'


fig, axes = plt.subplots(1, len(ratios), gridspec_kw={'width_ratios': ratios}, figsize=(12, 5), subplot_kw={'projection':'mantid'}) 


ax1 = axes[0]
ax1.pcolormesh(mtd[foldedSEDNames[0]], **colormesh_pars)
ax1.set_ylabel('E (meV)')
ax1.set_xlabel('')
ax1.tick_params(direction='in')

ax2 = axes[1]
ax2.pcolormesh(mtd[foldedSEDNames[1]], **colormesh_pars)
ax2.set_xlabel('')
ax2.get_yaxis().set_visible(False)
ax2.tick_params(direction='in')

ax3 = axes[2]
ax3.pcolormesh(mtd[foldedSEDNames[2]], **colormesh_pars)
ax3.set_xlabel('')
ax3.get_yaxis().set_visible(False)
ax3.tick_params(direction='in')

ax4 = axes[3]
ax4.pcolormesh(mtd[foldedSEDNames[3]], **colormesh_pars)
ax4.set_xlabel('')
ax4.get_yaxis().set_visible(False)
ax4.tick_params(direction='in')

axL = axes[-1]
mappable = axL.pcolormesh(mtd[foldedSEDNames[-1]], **colormesh_pars)
axL.get_yaxis().set_visible(False)
axL.set_xlabel('')
axL.tick_params(direction='in')


x0 = [ax.get_xlim() for ax in axes]

print(x0)

#ax1.set_xlim([x0[0], x0[0]])
ax1.set_xticks([x0[0][0], x0[0][1]])
ax1.set_xticklabels([r'$\Gamma$',''])

#ax2.set_xlim([x0[1], x0[1]])
ax2.set_xticks([x0[1][0], x0[1][1]])
#ax2.set_xticklabels([user_defined_path['path'][1][0],''])
ax2.set_xticklabels([user_defined_path['path'][1][0],''])

#ax3.set_xlim([x0[2], x0[2]])
ax3.set_xticks([x0[2][0], x0[2][1]])
ax3.set_xticklabels([user_defined_path['path'][2][0],''])

#ax4.set_xlim([x0[3], x0[3]])
ax4.set_xticks([x0[3][0], x0[3][1]])
ax4.set_xticklabels([user_defined_path['path'][3][0],''])

#axL.set_xlim([x0[4], x0[4]])
axL.set_xticks([x0[4][0], x0[4][1]])
#axL.set_xticklabels([user_defined_path['path'][4][0], user_defined_path['path'][4][1]])
axL.set_xticklabels([r'$\Gamma$', user_defined_path['path'][4][1]])

fig.suptitle('Folded Bi 300 K, 20 meV',  fontsize=22)

plt.subplots_adjust(wspace=0, hspace=0)

# Create the colorbar
cb_ax = fig.add_axes([.91,.124,.02,.754])
cbar = fig.colorbar(mappable,orientation='vertical',cax=cb_ax)
cbar.ax.tick_params(labelsize=15)

plt.rcParams.update({'font.size': 10})
fig.savefig(outputs_folder+'Bi_300K_20meV_SED.png')
