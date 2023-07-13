import os
import numpy as np

########################################################################################################
# Define list of dictionaries, each describing data (run) sets to be combined in a single mde workspace 
# Authors: A. Savici, I. Zaliznyak, March 2019.
########################################################################################################
def define_data_set(**kwargs):
    shared_folder='/SNS/ARCS/IPTS-21211/shared/Dipanshu/'
    raw_data_folder='/SNS/ARCS/IPTS-21211/nexus/'
    mde_folder='/SNS/ARCS/IPTS-21211/shared/Aiden/300K/mde/'

    data_set_list=[]
    # T=300K Ei=12meV
    data_set={'Runs':list(range(112488,112847+1)), #300K quick
              #'Runs':list(range(314383,314400+1)),    #List of runs, or list of lists of runs that are added together
              'BackgroundRuns':None,   #range(297325,297337)Options: None;list of runs that are added together
              'RawDataFolder':raw_data_folder,      #Options:raw_data_folder string
              'MdeFolder':mde_folder,               #Options:mde_folder string
              'MdeName':'FeSi_300K_Ei70meV_full',   #Options:mde_name string
              'BackgroundMdeName':'',      #Options:None;bg_mde_name string
              'MaskingDataFile':shared_folder+'van108467_new.nxs',         #Options:None;data_file_name
              'NormalizationDataFile':shared_folder+'van108467_new.nxs',   #Options:None;data_file_name
              'SampleLogVariables':{'OmegaMotorName':None,'Temperature':300.0},   #Options:None;LogVariableName;number
              'UBSetup':{'a':4.449,'b':4.449,'c':4.449,'alpha':90,'beta':90,'gamma':90,'u':'-0.0562,-0.0797,1','v':'0.9809,1,0.1348'}, #by hand
               #Data reduction options
              'Ei':70.,                             #Options: None;Ei_somehow_determined
              'T0': None,                           #Options: None;T0_determined_from_mantid
              'BadPulsesThreshold':None,            #Options: None;bg_pulses_threshold value
              'TimeIndepBackgroundWindow':None,     #Options: None;[Tib_min,Tib_max]
              'E_min':None,                         #Options: None;Emin
              'E_max':None,                         #Options: None;Emax
              'AdditionalDimensions':None,          #Options: None;list of triplets ("name", min, max)
              'EfCorrectionFunction':None,          #Options:None;'HYSPEC_default_correction';Custom_Ef_Correction_Function_Name
              }
    data_set_list.append(data_set)

    return data_set_list

