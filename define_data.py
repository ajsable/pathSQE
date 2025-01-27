import os
import numpy as np

########################################################################################################
# Define list of dictionaries, each describing data (run) sets to be combined in a single mde workspace 
# Authors: A. Savici, I. Zaliznyak, March 2019.
########################################################################################################
def define_data_set(**kwargs):
    shared_folder='/SNS/CNCS/IPTS-31965/shared/autoreduce/'
    raw_data_folder='/SNS/CNCS/IPTS-31965/nexus/'
    mde_folder='/SNS/CNCS/IPTS-31965/shared/chengjie/merged_mde/'

    T = 2
    E = 20

    data_set_list=[]
    data_set={'Runs':None, 
              'BackgroundRuns':None,   #range(297325,297337)Options: None;list of runs that are added together
              'RawDataFolder':raw_data_folder,      #Options:raw_data_folder string
              'MdeFolder':mde_folder,               #Options:mde_folder string
              'MdeName':'Bi_{}K_{}meV'.format(T,E), 
              'BackgroundMdeName':None,      #Options:None;bg_mde_name string
              'MaskingDataFile':shared_folder+'van_449014.nxs',         #Options:None;data_file_name; van_449014_unmasked_negative_detectors.nxs
              'NormalizationDataFile':shared_folder+'van_449014.nxs',   #Options:None;data_file_name
              'SampleLogVariables':{'OmegaMotorName':None,'Temperature':T},   #Options:None;LogVariableName;number FROM CHENGJIE
              'UBSetup':{'a':4.54,'b':4.54,'c':11.862,'alpha':90,'beta':90,'gamma':120,'u':'-1.14875,-0.12355,11.27797','v':'3.70658,0.06694,3.65633'}, #by hand FROM CHENGJIE
               #Data reduction options
              'Ei':None,                             #Options: None;Ei_somehow_determined
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

