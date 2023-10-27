import matplotlib.pyplot as plt
from mantid import plots
import sys
sys.path.append('/SNS/groups/dgs/DGS_SC_scripts')
from reduce_data_to_MDE import *
from slice_utils import *
import define_data
import define_slices
import os
import imp
from matplotlib.ticker import AutoLocator

directory_of_this_script = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, directory_of_this_script)

def set_int_ticks(int_ticks='False'):
    if int_ticks=='True':
        qdim=0
        for dim in range(3): # check to see how many and which q dim have 3 args, and set int ticks based on qdim bounds
            l = []
            for t in slice_description['Dimension{}Binning'.format(dim)].split(','):
                try:
                    l.append(float(t))
                except ValueError:
                    pass
            if len(l)==3 and qdim==0:
                plt.xticks(np.arange(np.floor(l[0]), np.ceil(l[2]+1), step=2))
                plt.xticks(np.arange(np.floor(l[0]), np.ceil(l[2]+1), step=1), minor=True)
                plt.xlim([l[0], l[2]])
                qdim=1
            elif len(l)==3 and qdim==1:
                plt.yticks(np.arange(np.floor(l[0]), np.ceil(l[2]+1), step=2))
                plt.yticks(np.arange(np.floor(l[0]), np.ceil(l[2]+1), step=1), minor=True)
                plt.ylim([l[0], l[2]])
        return
    else:
        return
        
#############################################

T_vec = [650]#[10,300,650]
E_vec = [70]#[70]

#############################################
# MAIN PROGRAM
#############################################
imp.reload(define_data)
imp.reload(define_slices)

datasets=define_data.define_data_set(T_vec,E_vec)

for i in range(len(datasets)):
    ds=datasets[i]
    reduce_data_to_MDE([ds])
    
    T_curr = T_vec[i//len(E_vec)]
    E_curr = E_vec[(i)%len(E_vec)]
    slice_descriptions=define_slices.make_slices_aiden(T_curr,E_curr)
    filename = directory_of_this_script + '/plots/'

    for slice_description in slice_descriptions:
        #plt.rcParams.update({'font.size': 10})
        make_slice(ds,slice_description, solid_angle_ws=None, ASCII_slice_folder='')
        fig,ax=plt.subplots(subplot_kw={'projection':'mantid'})
        plt.rcParams['figure.dpi'] = 300
        plt.rcParams['savefig.dpi'] = 300
        plt.set_cmap('viridis')
        c=plot_slice(slice_description, ax=ax, cbar_label='Intensity')
        set_axes_parameters(ax,**slice_description['Axes_parameters'])
        #ax.set_yscale('log')
        #ax.grid(True, which = 'both', axis = 'both')
       
        # make integer tick marks for non-energy axes for readability
        set_int_ticks(int_ticks=slice_description['Int_ticks'])
        
        # skip plots that are essentially empty
        if np.sum(mtd[slice_description['Name']].getNumEventsArray()) > 100:
            #fig.savefig(filename+slice_description['Name']+'.png')
            #fig.show()
        else:
            print("The plot is empty.")
        
