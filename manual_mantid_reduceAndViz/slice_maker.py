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
directory_of_this_script = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, directory_of_this_script)

#############################################
# MAIN PROGRAM
#############################################
imp.reload(define_data)
imp.reload(define_slices)

datasets=define_data.define_data_set()
ds=datasets[0]
reduce_data_to_MDE([ds])

#SetUB(ds['MdeName'].strip(), a=3.95,b=3.95,c=7.6,alpha=90,beta=90,gamma=90,u='1,1,0',v='0,0,1')
#SetUB(ds['MdeName'].strip(), a=3.95,b=3.95,c=7.6,alpha=90,beta=90,gamma=90,u='-0.41307364752,0.369818133844,-0.832228760381',v='-0.657979223253,0.510606268289,0.553484038209')
#SetUB(ds['MdeName'].strip(), a=3.95,b=3.95,c=7.6,alpha=90,beta=90,gamma=90,u='0.170384655046,-0.333251568204,-0.927314650814',v='0.269178063905,0.921022832833,-0.281531723456')
#SetUB(ds['MdeName'].strip(), a=3.95,b=3.95,c=7.6,alpha=90,beta=90,gamma=90,u='-0.336, -0.336, 0.94157',v='0.94157,0.94157,-0.3368')
#SetUB(ds['MdeName'].strip(), a=3.95,b=3.95,c=7.6,alpha=90,beta=90,gamma=90,u='0.42, 0.42, 0.91',v='0.91,0.91,0.42')
#SetUB(ds['MdeName'].strip(), a=3.95,b=3.95,c=7.6,alpha=90,beta=90,gamma=90,u='0.431, 0.431, 0.903',v='0.903,0.903,0.431') #by hand rotating, looks pretty good
#SetUB(ds['MdeName'].strip(), a=3.95,b=3.95,c=7.6,alpha=90,beta=90,gamma=90,u='0.402707790974,0.443662336676,0.800618614638',v='0.222448996793,-0.895892946423,0.384567643432')
#SetUB(ds['MdeName'].strip(), a=4.02,b=3.85,c=7.6,alpha=90,beta=90,gamma=90,u='1,1,0',v='0,0,1')
#SetUB(ds['MdeName'].strip(), a=3.9,b=3.9,c=7.6,alpha=90,beta=90,gamma=90,u='-0.6869, -0.7199, 0.09954',v='-0.0277, 0.0628, 0.99764 ')
#SetUB(ds['MdeName'].strip(), a=3.9,b=3.9,c=7.6,alpha=90,beta=90,gamma=90,v='-0.6869, -0.7199, 0.09954',u='-0.0277, 0.0628, 0.99764 ')


#slice_descriptions=define_slices.map_3D_loop()
#slice_descriptions=define_slices.slices_loop3()
#slice_descriptions=define_slices.slices_perp_rod2()
#slice_descriptions=define_slices.maps_0p5()
#slice_descriptions=define_slices.maps()
#slice_descriptions=define_slices.slices_0k0()
slice_descriptions=define_slices.make_slices_aiden()



mde_data=define_data.define_data_set()
reduce_data_to_MDE(mde_data)

#print (directory_of_this_script)
filename = directory_of_this_script


#print(directory_of_this_script)

#slice_description=slice_descriptions[0] #take only one slice
#for slice_description in slice_descriptions[-1:]:
for slice_description in slice_descriptions:
    make_slice(ds,slice_description, solid_angle_ws=None, ASCII_slice_folder='')
    #fig,ax=plt.subplots(subplot_kw={'projection':'mantid'})
    ## MY EDITS ##
    #plt.rcParams['figure.dpi'] = 300
    #plt.rcParams['savefig.dpi'] = 300
    #plt.set_cmap('jet')
    ##############
    #c=plot_slice(slice_description, ax=ax, cbar_label='Intensity')
    #set_axes_parameters(ax,**slice_description['Axes_parameters'])
    #fig.savefig(filename+slice_description['Name']+'.png')

    #fig.show()    
#for slice_description in slice_descriptions[-1:]:
#for slice_description in slice_descriptions:
#    make_slice(ds,slice_description, solid_angle_ws=None, ASCII_slice_folder=filename)