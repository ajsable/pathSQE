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



outputs_folder = '/SNS/ARCS/IPTS-13861/shared/Aiden/comprehensive/pathSQE_testing_FeSi/10K_40meV_foldWSims_paper/'

prim2mantid_transMatrix = np.array([[1,0,0], [0,1,0], [0,0,1]])
user_defined_path={'path':[ ('Gamma','X'), ('X','M'), ('M','Gamma'), ('Gamma','R'), ('R','X') ],  
                    '1d_points':['Gamma', 'R', 'X'],
                    'point_coords':{ 'Gamma': [0.0, 0.0, 0.0], 'X': [0.0, 0.5, 0.0], 'M': [0.5, 0.5, 0.0], 'R': [0.5, 0.5, 0.5] }}        
T = 10
Ebins = '0,0.25,40'


data_path = outputs_folder+'folding_progess/folded_sim_path_plot_96.npz'
with np.load(data_path, allow_pickle=True) as data_npz:
    data_arrays = [data_npz[key].T for key in data_npz]  # Extract all arrays


ratios = []
for seg in user_defined_path['path']:    
    pt1_prim = np.array(user_defined_path['point_coords'][seg[0]])
    pt2_prim = np.array(user_defined_path['point_coords'][seg[1]])
    norm_prim = np.linalg.norm(pt2_prim-pt1_prim)

    ratios.append(norm_prim)



plt.rcParams['figure.dpi'] = 600
plt.rcParams['savefig.dpi'] = 600
plt.rcParams.update({'font.size': 28})

colormesh_pars={}
colormesh_pars['norm']=SymLogNorm(linthresh=1e-1,vmin=1e-1,vmax=1e2) # data (linthresh=1e-3,vmin=1e-3,vmax=2e-1)
colormesh_pars['cmap']='viridis'


fig, axes = plt.subplots(1, len(ratios), gridspec_kw={'width_ratios': ratios}, figsize=(12, 5), subplot_kw={'projection':'mantid'}) 


ax1 = axes[0]
ax1.pcolormesh(data_arrays[0], **colormesh_pars)
ax1.set_ylabel('E (meV)')
ax1.set_xlabel('')
ax1.tick_params(direction='in')

ax2 = axes[1]
ax2.pcolormesh(data_arrays[1], **colormesh_pars)
ax2.set_xlabel('')
ax2.get_yaxis().set_visible(False)
ax2.tick_params(direction='in')

ax3 = axes[2]
ax3.pcolormesh(data_arrays[2], **colormesh_pars)
ax3.set_xlabel('')
ax3.get_yaxis().set_visible(False)
ax3.tick_params(direction='in')

ax4 = axes[3]
ax4.pcolormesh(data_arrays[3], **colormesh_pars)
ax4.set_xlabel('')
ax4.get_yaxis().set_visible(False)
ax4.tick_params(direction='in')

axL = axes[-1]
mappable = axL.pcolormesh(data_arrays[-1], **colormesh_pars)
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
ax4.set_xticklabels([r'$\Gamma$',''])

#axL.set_xlim([x0[4], x0[4]])
axL.set_xticks([x0[4][0], x0[4][1]])
#axL.set_xticklabels([user_defined_path['path'][4][0], user_defined_path['path'][4][1]])
axL.set_xticklabels([user_defined_path['path'][4][0], user_defined_path['path'][4][1]])

fig.suptitle('FeSi sim 100 K, 40 meV',  fontsize=22)

plt.subplots_adjust(wspace=0, hspace=0)

# Create the colorbar
cb_ax = fig.add_axes([.91,.124,.02,.754])
cbar = fig.colorbar(mappable,orientation='vertical',cax=cb_ax)
cbar.ax.tick_params(labelsize=15)

plt.rcParams.update({'font.size': 10})
fig.savefig(outputs_folder+'FeSi_10K_40meV_sim.png')
