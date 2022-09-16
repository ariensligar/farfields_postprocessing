from calendar import EPOCH
import pickle
import sys
import os

import sys
import pathlib

from matplotlib import pyplot as plt, cm

root_dir1 = pathlib.Path('__file__').absolute().parent
sys.path.append(str(root_dir1)) #path where utilitiesHfss.py script sits

from utilitiesHfss import Load_FF_Fields
import numpy as np
import math
import logging


log_dir = './log/'
if not os.path.exists(log_dir):
    os.makedirs(log_dir)
logging.basicConfig(filename=f'{log_dir}example88.log', level=logging.DEBUG)

global PyVectorDrivingAntennaGainC17_init
global PyVectorDrivingAntennaGainC17_Inputs
global PyVectorDrivingAntennaGainC17_Outputs
global all_eep
global lattice_vectors
global all_qtys
global oldPhi
global oldTheta
global gen_image
global oldEpoch
global counter
oldPhi = 0
oldTheta = 0
gen_image = False
oldEpoch = 0
counter = 1
all_eep = None

PyVectorDrivingAntennaGainC17_init = -1


def is_single_phase_reference(location_dict):
    empty_vec = np.zeros((len(location_dict), 3))
    for n, loc in enumerate(location_dict):
        test = location_dict[loc]
        empty_vec[n] = test
    if len(set(list(empty_vec[:, 0]))) == 1 and len(set(list(empty_vec[:, 1]))) == 1 and len(
            set(list(empty_vec[:, 2]))) == 1:
        return True
    else:
        return False

#==========================================================================
# PY_CalcObject() fctn
#==========================================================================
def PyVectorDrivingAntennaGainC17 ( argList ):
    logging.debug('arguments list')
    callMode = str(argList[0])
    logging.debug(f'Call Mode: {callMode}')
    if callMode == 'None':

        retVal = PyVectorDrivingAntennaGainC17_compute( argList )    # do compute
    elif callMode == 'register' :
        global PyVectorDrivingAntennaGainC17_init
        PyVectorDrivingAntennaGainC17_init = -1
        retVal = PyVectorDrivingAntennaGainC17_register()
    elif callMode == 'compute' :
        retVal = PyVectorDrivingAntennaGainC17_compute( argList )    # do compute
    else:
        retVal = []    # # bad call, return empty list
    return retVal

def PyVectorDrivingAntennaGainC17_register():
    logging.debug('Registering Script')
    return [
        ["ArgumentType = Output", "Name = AntennaGain", "ArgumentName = AntennaGain"],
        ["ArgumentType = Output", "Name = Beamwidth", "ArgumentName = Beamwidth"],
        ["ArgumentType = Output", "Name = AntennaMaxGain", "ArgumentName = AntennaMaxGain"],
        ["ArgumentType = Output", "Name = IntegratedGain", "ArgumentName = IntegratedGain"],
        ["ArgumentType = Output", "Name = AntennaCoordSystem", "ArgumentName = AntennaCoordSystem"],
        ["ArgumentType = Output", "Name = DynamicGain", "ArgumentName = Value"],
        ["ArgumentType = Input", "Name = DateUTC", "ArgumentName = DateUTC", "Type = Value"],
        ["ArgumentType = Input", "Name = CbName", "ArgumentName = CbName", "Type = Value"],
        ["ArgumentType = Input", "Name = Frequency", "ArgumentName = Frequency", "Type = Value"],
        ["ArgumentType = Input", "Name = AzimuthAngle", "ArgumentName = AzimuthAngle", "Type = Value"],
        ["ArgumentType = Input", "Name = ElevationAngle", "ArgumentName = ElevationAngle", "Type = Value"],
        ["ArgumentType = Input", "Name = AntennaPosLLA", "ArgumentName = AntennaPosLLA", "Type = Value"],
        ["ArgumentType = Input", "Name = AntennaCoordSystem", "ArgumentName = AntennaCoordSystem", "Type = Value"],
        ["ArgumentType = Input", "Name = PhiPointing", "ArgumentName = PhiPointing", "Type = Angle"],
        ["ArgumentType = Input", "Name = ThetaPointing", "ArgumentName = ThetaPointing", "Type = Angle"]
    ]



def PyVectorDrivingAntennaGainC17_compute( inputData ):
    # NOTE: argList[0] is the call Mode, which is either None or 'compute'
    global PyVectorDrivingAntennaGainC17_init
    global PyVectorDrivingAntennaGainC17_Inputs
    global all_eep
    global lattice_vectors
    global all_qtys
    global oldPhi
    global oldTheta
    global phi_scan
    global theta_scan
    global gen_image
    global oldEpoch
    global hfss_utils
    global counter

    logging.debug(f'Number of times called, {counter}')
    counter += 1
    ##############################################################################
    #
    #   USER INPUTS
    #
    ##############################################################################
    #load results from file, instead of loading them dynamiclaly from HFSS
    eep_results_file = r'C:\Users\asligar\OneDrive - ANSYS, ' \
                    r'Inc\Documents\Scripting\github\farfields_postprocessing\example_projects\ffd_data\exportelement.txt '
    lattice_vectors = None
    ##############################################################################
    #
    #   END USER INPUTS
    #
    ##############################################################################

    script_path = os.getcwd()
    full_path = os.path.abspath(eep_results_file)
    ffd_base_path, file_name = os.path.split(full_path)
    split_ext_filename = os.path.splitext(eep_results_file)

    if PyVectorDrivingAntennaGainC17_init != 1:
        if split_ext_filename[1] is 'eep':
            f = open(eep_results_file, 'rb')
            obj = pickle.load(f)
            f.close()
            locals().update(obj)
            all_eep = obj['all_eep']
            lattice_vectors = obj['lattice_vectors']
            logging.debug('Loading far field data from eep file')
            hfss_utils = Load_FF_Fields(all_eep, lattice_vectors=lattice_vectors)
            logging.debug('EEP file Loaded')
            logging.debug(f'File Type: EEP, {eep_results_file}')

        else:
            logging.debug(f'File Type: FFD, {eep_results_file}')
            with open(eep_results_file, 'r') as reader:
                lines = [line.split(None) for line in reader]
            reader.close()
            lines = lines[1:]  # remove header

            path_dict = {}
            location_dict = {}
            for pattern in lines:
                if len(pattern) >= 2:
                    port = pattern[0]
                    if ':' in port:
                        port = port.split(':')[0]
                    path_dict[port] = ffd_base_path + '/' + pattern[1] + '.ffd'
                if len(pattern) == 5:  # it contains position information
                    x = float(pattern[2])
                    y = float(pattern[3])
                    z = float(pattern[4])
                    xyz = np.array([x, y, z])
                    location_dict[port] = xyz

            is_single_phase_reference_used = is_single_phase_reference(location_dict)
            if is_single_phase_reference_used:
                lattice_vectors = None
                logging.debug('FFD files all referenced to the same location, please provide lattice vectors')
                return
            else:
                lattice_vectors = None
                logging.debug('FFD contain reference location, using these positions for array processing')

            logging.debug('Loading far field data from ffd files')
            hfss_utils = Load_FF_Fields(path_dict, lattice_vectors=None, element_position_dict=location_dict)
            logging.debug(f'FFD Files Loaded')
            logging.debug(f'File Type: FFD, {eep_results_file}')


    logging.debug('Computing')

    if PyVectorDrivingAntennaGainC17_init < 0:
        logging.debug(f'Initializing plug in')
        #log.write(time.ctime() + " - Initializing plugin\n")
        phi_scan = -9999
        theta_scan = -9999
        freq = -9999

        PyVectorDrivingAntennaGainC17_init = 1
        PyVectorDrivingAntennaGainC17_Inputs = g_PluginArrayInterfaceHash['PyVectorDrivingAntennaGainC17_Inputs']
        PyVectorDrivingAntennaGainC17_Outputs = g_PluginArrayInterfaceHash['PyVectorDrivingAntennaGainC17_Outputs']

    phiPointingRad = float(inputData[PyVectorDrivingAntennaGainC17_Inputs['PhiPointing']])
    thetaPointingRad = float(inputData[PyVectorDrivingAntennaGainC17_Inputs['ThetaPointing']])
    az = float(inputData[PyVectorDrivingAntennaGainC17_Inputs['AzimuthAngle']])
    el = float(inputData[PyVectorDrivingAntennaGainC17_Inputs['ElevationAngle']])
    freq = float(inputData[PyVectorDrivingAntennaGainC17_Inputs['Frequency']])
    evalTime = inputData[PyVectorDrivingAntennaGainC17_Inputs['DateUTC']]
    #############################################################################################
    # USER ANTENNA GAIN MODEL AREA.
    # PLEASE REPLACE THE CODE BELOW WITH YOUR ANTENNA GAIN COMPUTATION MODEL
    #############################################################################################
    # NOTE: the outputs that are returned MUST be in the same order as registered
    # AntennaGain (dB), gain of the antenna at time and in the Azi-Elev direction off the boresight.
    # Beamwidth (Rad) is the 3-dB beamwith of the antenna.
    # AntennaMaxGain (dB) is the maximum ( possibly boresight gain of the antenna)
    # IntegratedGain of the antenna (range 0-1) used for antenna Noise computation.
    #
    """""""""""""""""""""
    START OF ANSYS SCRIPT
    """""""""""""""""""""



    phi_scan = phiPointingRad*(180/math.pi)
    theta_scan = thetaPointingRad*(180/math.pi)


    hfss_utils.freq = freq
    if phi_scan != oldPhi and theta_scan != oldTheta:
        #log.write("------------- Getting gain pattern for phi: " + str(phi_scan) + " theta: " + str(theta_scan) +"\n")
        logging.debug('Beamforming Started')
        all_qtys = hfss_utils.beamform(
            phi_scan=phi_scan,
            theta_scan=theta_scan)
        logging.debug('Beamforming Completed')
        oldPhi = phi_scan
        oldTheta = theta_scan

    azDeg = float(az)*(180/math.pi)
    elDeg = float(el)*(180/math.pi)
    azIndx = round(azDeg)
    elIndx = round(elDeg)
    th = all_qtys['Theta']
    ph = all_qtys['Phi']

    phiIndex = min(range(len(ph)), key=lambda i: abs(ph[i]-azDeg))
    thetaIndex = min(range(len(th)), key=lambda i: abs(th[i]-elDeg))
    eff = 0.55
    dia = 1.0
    lambda_value = 299792458.0 / freq
    thetab = lambda_value / (dia * math.sqrt(eff))
    gainArray = 10*np.log10(np.abs(all_qtys['RealizedGain']))

    gain = float(gainArray[thetaIndex, phiIndex])
    gmax =np.max(gainArray)
    logging.debug(f'Max Gain {gmax}')
    temp_index = np.argmax(gainArray)
    temp_theta = all_qtys['Theta'][np.argmax(all_qtys['RealizedGain'])%91]
    temp_phi = all_qtys['Phi'][np.argmax(all_qtys['RealizedGain'])%181]
    logging.debug(f'Theta: {temp_theta}')
    logging.debug(f'Phi: {temp_phi}')

    logging.info(f'{oldEpoch} - {evalTime}')
    logging.info(oldEpoch != evalTime)
    if oldEpoch != evalTime:
        logging.info(f'{oldEpoch} - {evalTime}')
        gen_image = False

    if gen_image:

        save_name = f'3D_Polar_Plot_{temp_theta}_{temp_phi}'
        logging.debug(f'Saving Plot of Far Field Patter: {save_name}')
        polar_plot_3d(all_qtys,
                       save_name =save_name,
                       save_plot=False,
                       show_plot=False,
                       output_path = './plots/',
                       folder_path = './plots/',
                       dB=True,
                       multiple_angles = False)
    oldEpoch = evalTime
    gen_image = False

    """""""""""""""""""""
    END OF ANSYS SCRIPT
    """""""""""""""""""""
    return [
        gain,                   # PyVectorDrivingAntennaGainC17_Outputs.AntennaGain (returns gain based on azimuth and elevation)
        thetab,                    # PyVectorDrivingAntennaGainC17_Outputs.Beamwidth
        gmax,                   # PyVectorDrivingAntennaGainC17_Outputs.AntennaMaxGain (taken from Ansys file)
        .5,                    # PyVectorDrivingAntennaGainC17_Outputs.IntegratedGain
        0,                        # PyVectorDrivingAntennaGainC17_Outputs.AntennaCoordSystem, AntennaCoordSystem return 0 for Polar and 1 for Rectangular
        1,                       #dynamic output
    ]

def polar_plot_3d(data,
                    save_name ="3D_Polar_Plot_Envelope",
                    save_plot=True,
                    show_plot=True,
                    output_path = '',
                    folder_path = './',
                    dB=True,
                    multiple_angles = True):
    if dB:
        ff_data = 10*np.log10(data['RealizedGain'])
        #renormalize to 0 and 1
        ff_max_dB = np.max(ff_data)
        ff_min_dB = np.min(ff_data)
        ff_data_renorm = (ff_data-ff_min_dB)/(ff_max_dB-ff_min_dB)
    else:
        ff_data = data['RealizedGain']
        #renormalize to 0 and 1
        ff_max = np.max(ff_data)
        ff_min = np.min(ff_data)
        ff_data_renorm = (ff_data-ff_max)/(ff_max-ff_min)
    legend = []
    theta = np.deg2rad(np.array(data['Theta']))
    phi = np.deg2rad(np.array(data['Phi']))
    phi_grid,theta_grid = np.meshgrid(phi, theta)
    r = np.reshape(ff_data_renorm,(len(data['Theta']),len(data['Phi'])))
    x = r * np.sin(theta_grid) * np.cos(phi_grid)
    y = r * np.sin(theta_grid) * np.sin(phi_grid)
    z = r * np.cos(theta_grid)
    fig1 = plt.figure()
    ax1 = fig1.add_subplot(1, 1, 1, projection="3d")
    my_col = cm.jet(r/np.amax(r))
    plot = ax1.plot_surface(
        x, y, z, rstride=1, cstride=1, cmap=plt.get_cmap("jet"),facecolors = my_col, linewidth=0, antialiased=True, alpha=0.9)
    #fig1.set_size_inches(22.5, 22.5)
    plt.colorbar(plot)
    if save_plot:
       if multiple_angles:
           list_of_observations= [(0,0),(0,90),(0,180),(0,270),(90,0),(45,45),(45,-45),(-45,-45)]
           for n, observe in enumerate(list_of_observations):
               ax1.view_init(elev=observe[0], azim=observe[1])
               save_name = save_name + '_' + str(n) + '.png'
               save_name_full = folder_path + save_name
               plt.savefig(save_name_full,dpi=300)
               logging.info(f'Save Image: {save_name_full}')
       else:
            save_name_full = folder_path + save_name + '.png'
            plt.savefig(save_name_full,dpi=300)
    if show_plot:
        plt.show()

def main():
    if __name__ == "__main__":
        main()
