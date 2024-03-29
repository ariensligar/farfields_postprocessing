# -*- coding: utf-8 -*-
"""

Contains function for extracting far field values and embedded element patterns

Arien Sligar (arien.sligar@ansys.com)
Last Update: 06/24/2020
"""
import numpy as np
import os
import math
import pyvista as pv
import time as walltime
import pickle

import logging
log = logging.getLogger(__name__)






def get_array_index(port_name):

    str1 = port_name.split('[', 1)[1].split(']', 1)[0]
    index_list = str1.split(',')
    index_list = [int(i) for i in index_list]
    #logging.debug(f'Array Index: {index_list}')
    return index_list


def calc_relative_phase(pos, wl, theta, phi):
    phaseConstant = (2 * np.pi / wl)

    xVector = -pos[0] * np.sin(theta) * np.cos(phi)
    yVector = -pos[1] * np.sin(theta) * np.sin(phi)
    zVector = -pos[2] * np.cos(theta)

    phaseOfIncidentWaveAtElement = phaseConstant * (xVector + yVector + zVector)

    return phaseOfIncidentWaveAtElement


def offset_antenna_phase(re_theta, re_phi, pos, freq, theta_range, phi_range):
    """
    beamforming is done based on a single reference point. If the array has been exported with a xyz offset
    this will transform the ffd file to the same reference position. This does not impact anything in the simulation
    only the visualization of the far fields. I could have also instead done the beamforming based on absolute position
    but in the end it is the same

    Parameters
    ----------
    re_theta : TYPE
        DESCRIPTION.
    re_phi : TYPE
        DESCRIPTION.
    pos : TYPE
        DESCRIPTION.
    freq : TYPE
        DESCRIPTION.
    theta_range : TYPE
        DESCRIPTION.
    phi_range : TYPE
        DESCRIPTION.

    Returns
    -------
    re_theta_updated : TYPE
        DESCRIPTION.
    re_phi_updated : TYPE
        DESCRIPTION.

    """
    logging.debug('Offsetting Phase To Common Reference')
    re_theta_updated = np.ones(re_theta.shape, dtype='complex')
    re_phi_updated = np.ones(re_phi.shape, dtype='complex')
    wl = 3e8 / freq
    idx = 0

    # k = (2 * np.pi / wl)
    # theta_grid, phi_grid = np.meshgrid(theta_range, phi_range, indexing='ij')
    # xv = -pos[0] * np.sin(theta_grid) * np.cos(phi_grid)
    # yv = -pos[1] * np.sin(theta_grid) * np.sin(phi_grid)
    # zv = -pos[2] * np.cos(theta_grid)
    # relative_phase = k * (xv + yv + zv)

    for theta_idx, theta in enumerate(theta_range):
        for phi_idx, phi in enumerate(phi_range):
            relativePhase = calc_relative_phase(pos, wl, np.deg2rad(theta), np.deg2rad(phi))
            re_theta_updated[idx] = np.abs(re_theta[idx]) * np.exp((np.angle(re_theta[idx]) - relativePhase) * 1j)
            re_phi_updated[idx] = np.abs(re_phi[idx]) * np.exp((np.angle(re_phi[idx]) - relativePhase) * 1j)
            idx += 1
    logging.debug('Offsetting Phase To Common Reference: Done')
    return re_theta_updated, re_phi_updated


class Load_FF_Fields():
    """
    Load and compute near field data as defined in the codebook through superpostion
    of the individual near field values for each port. Scaled by mag/phase
    """

    def __init__(self, ffd_dict, lattice_vectors=None, element_position_dict=None):

        time_before = walltime.time()
        print('Loading Embedded Element Patterns...')
        logging.debug('Loading Embedded Element Patterns..')
        self.data_dict = {}
        self.data_for_eep_export = {}
        self.freq = 1e9
        self.taper = 'flat'
        self.element_position_dict = element_position_dict
        self.position_by_index = {}
        valid_ffd = True
        all_ports = list(ffd_dict.keys())
        if os.path.exists(ffd_dict[all_ports[0]]):
            logging.debug(f'Loading: {ffd_dict[all_ports[0]]}')
            with open(ffd_dict[all_ports[0]], 'r') as reader:
                theta = [int(i) for i in reader.readline().split()]
                phi = [int(i) for i in reader.readline().split()]
                num_freq = int(reader.readline().split()[1])
                frequency = float(reader.readline().split()[1])
            reader.close()
            port_idx_all = []
            results_dict = {}
            for port in ffd_dict.keys():
                if ':' in port:
                    port = port.split(':')[0]
                temp_dict = {}
                theta_range = np.linspace(*theta)
                phi_range = np.linspace(*phi)

                ntheta = len(theta_range)
                nphi = len(phi_range)

                if os.path.exists(ffd_dict[port]):
                    port_idx = get_array_index(port)
                    port_idx_all.append(port_idx)
                    # self.position_by_index[port_idx] = self.element_position_dict[port]
                    eep_txt = np.loadtxt(ffd_dict[port], skiprows=4)
                    Etheta = np.vectorize(complex)(eep_txt[:, 0], eep_txt[:, 1])
                    Ephi = np.vectorize(complex)(eep_txt[:, 2], eep_txt[:, 3])

                    # eep=np.column_stack((etheta, ephi))
                    temp_dict['Theta'] = theta_range
                    temp_dict['Phi'] = phi_range
                    if lattice_vectors is not None:
                        temp_dict['rETheta'] = Etheta
                        temp_dict['rEPhi'] = Ephi
                    else:
                        # translate fields by actual position so they are placed at common reference position
                        # beamforming will later be based on this assumption
                        temp_dict['rETheta'], temp_dict['rEPhi'] = offset_antenna_phase(
                            Etheta,
                            Ephi,
                            element_position_dict[port],
                            frequency,
                            theta_range,
                            phi_range)

                    self.data_dict[port] = temp_dict
                    self.data_for_eep_export[port] = [theta_range, phi_range, temp_dict['rETheta'], temp_dict['rEPhi']]
                else:
                    valid_ffd = False

            # differential area of sphere, based on observation angle
            self.d_theta = np.abs(theta_range[1] - theta_range[0])
            self.d_phi = np.abs(phi_range[1] - phi_range[0])
            self.diff_area = np.radians(self.d_theta) * np.radians(self.d_phi) * np.sin(np.radians(theta_range))
            self.num_samples = len(temp_dict['rETheta'])
            self.all_port_names = list(self.data_dict.keys())
            self.solution_type = 'DrivenModal'
            if frequency:
                self.freq = frequency
                print(f"INFO: Using frequency from ffd file: {self.freq}")
                logging.debug(f'Using frequency from ffd file: {self.freq}')
        else:
            valid_ffd = False
            print('ERROR: Far Field Files are Missing')
            logging.debug('ERROR: Far Field Files are Missing')

        self.valid_ffd = valid_ffd
        self.lattice_vectors = lattice_vectors
        # if self.lattice_vectors is None:
        # array_idx = np.array(port_idx_all)
        # row_min = np.min(array_idx[:, 0])
        # row_max = np.max(array_idx[:, 0])
        # col_min = np.min(array_idx[:, 1])
        # col_max = np.max(array_idx[:, 1])

        # a = np.empty((row_max,col_max,))

        # lattice_vectors = get_lattice_vectors_from_positions(self.element_position_dict)
        if self.lattice_vectors is not None:
            self.Ax = float(lattice_vectors[0])
            self.Ay = float(lattice_vectors[1])
            self.Bx = float(lattice_vectors[3])
            self.By = float(lattice_vectors[4])
            logging.debug(f'Lattice Vectors: {lattice_vectors}')

        elapsed_time = walltime.time() - time_before
        print(f'Loading Embedded Element Patterns...Done: {elapsed_time}seconds')
        logging.debug(f'Loading Embedded Element Patterns...Done: {elapsed_time}seconds')

    def find_array_min_max_idx(self, all_port_names=None):
        """
        Parameters
        ----------
        all_port_names : list of strings, if this isn't suplied we will assume
        it was already intialized in the main script before calling ff_beamsteer()
            DESCRIPTION. The default is None.

        Returns
        -------
        [int, int, int, int]
            [column min, column max, row min, row max] of array

        """
        if all_port_names is None:
            all_port_names = self.all_port_names
        row_min = 1
        row_max = 1
        col_min = 1
        col_max = 1
        rows = []
        cols = []
        for port_string in all_port_names:
            index_str = get_array_index(port_string)
            rows.append(int(index_str[1]))
            cols.append(int(index_str[0]))

        row_min = np.min(rows)
        row_max = np.max(rows)
        col_min = np.min(cols)
        col_max = np.max(cols)
        return [col_min, col_max, row_min, row_max]

    def FindArrayCenterAndEdge(self):
        '''
        Find teh center and edge of our array, assumes all ports in far field
        mapping file are active ports. Call before AssignWeight()

        Returns
        -------
        None. Sets Rmax, Amax... etc

        '''
        AMax = 0
        BMax = 0
        RMax = 0
        XMax = 0
        YMax = 0
        CenterA = 0
        CenterB = 0
        CenterX = 0
        CenterY = 0

        # collecting all active cells inside the specified region
        activeCells = []

        for i in range(0, len(self.all_port_names)):
            index_str = get_array_index(self.all_port_names[i])
            row = int(index_str[1])
            col = int(index_str[0])
            a = row
            b = col

            activeCells.append((a, b))  # because ffd is assuming all ffd files are active
        if len(activeCells) == 0: return

        [a_min, a_max, b_min, b_max] = self.find_array_min_max_idx()
        CenterA = (a_min + a_max) / 2
        CenterB = (b_min + b_max) / 2

        if self.lattice_vectors is None:
            all_pos = []
            for port in self.element_position_dict:
                all_pos.append(self.element_position_dict[port])
            all_pos = np.array(all_pos)
            x_min = np.min(all_pos[:, 0])
            x_max = np.max(all_pos[:, 0])
            y_min = np.min(all_pos[:, 1])
            y_max = np.max(all_pos[:, 1])

            CenterX = (x_max - x_min) / 2
            CenterY = (y_max - y_min) / 2
        else:
            CenterX = (CenterA + 0.5) * self.Ax + (CenterB + 0.5) * self.Bx
            CenterY = (CenterA + 0.5) * self.Ay + (CenterB + 0.5) * self.By

        self.CenterA = CenterA
        self.CenterB = CenterB
        self.CenterX = CenterX
        self.CenterY = CenterY
        # find the distance from the edge to the center
        AMax = a_max - a_min
        BMax = b_max - b_min

        self.AMax = AMax
        self.BMax = BMax

        if self.lattice_vectors is None:
            for port in self.element_position_dict:
                x = self.element_position_dict[port][0]
                y = self.element_position_dict[port][1]
                x_dis = abs(x - CenterX)
                y_dis = abs(y - CenterY)
                distance = math.sqrt(x_dis ** 2 + y_dis ** 2)
                XMax = max(XMax, x_dis)
                YMax = max(YMax, y_dis)
                RMax = max(RMax, distance)
        else:
            for a, b in activeCells:
                x = (a + 0.5) * self.Ax + (b + 0.5) * self.Bx
                y = (a + 0.5) * self.Ay + (b + 0.5) * self.By
                x_dis = abs(x - CenterX)
                y_dis = abs(y - CenterY)
                distance = math.sqrt(x_dis ** 2 + y_dis ** 2)
                XMax = max(XMax, x_dis)
                YMax = max(YMax, y_dis)
                RMax = max(RMax, distance)

        self.RMax = RMax
        self.XMax = XMax
        self.YMax = YMax
        self.RMax *= 2
        self.XMax *= 2
        self.YMax *= 2

    def ElementPosition(self, a, b):
        a = int(a)
        b = int(b)

        x = (a + 0.5) * self.Ax + (b + 0.5) * self.Bx
        y = (a + 0.5) * self.Ay + (b + 0.5) * self.By
        x_dis = x - self.CenterX
        y_dis = y - self.CenterY

        return np.array([x_dis, y_dis, 0])

    def AssignWeight(self, a, b, taper='flat'):
        '''
        assign weights to our array

        Parameters
        ----------
        a : INT
            index of array, column.
        b : INT
            index or array, row.
        taper : string, optional
            DESCRIPTION. This is teh type of taper we want to apply. The default is 'flat'.

        Returns
        -------
        float
            weight to applied to specific index of array.

        '''
        a = int(a)
        b = int(b)
        if taper.lower() == 'flat':  # Flat
            return 1

        cosinePow = 1
        edgeTaper_dB = -200

        edgeTaper = 10 ** ((float(edgeTaper_dB)) / 20)

        threshold = 1e-10
        length_in_direction1 = 0
        max_length_in_dir1 = 0
        length_in_direction2 = 0
        max_length_in_dir2 = 0
        w1 = w2 = None

        # find the distance between current cell and array center in terms of index
        length_in_direction1 = a - self.CenterA
        length_in_direction2 = b - self.CenterB
        max_length_in_dir1 = self.AMax
        max_length_in_dir2 = self.BMax

        if taper.lower() == 'cosine':  # Cosine
            if max_length_in_dir1 < threshold:
                w1 = 1
            else:
                w1 = (1 - edgeTaper) * (
                    math.cos(math.pi * length_in_direction1 / max_length_in_dir1)) ** cosinePow + edgeTaper
            if max_length_in_dir2 < threshold:
                w2 = 1
            else:
                w2 = (1 - edgeTaper) * (
                    math.cos(math.pi * length_in_direction2 / max_length_in_dir2)) ** cosinePow + edgeTaper
        elif taper.lower() == 'triangular':  # Triangular
            if max_length_in_dir1 < threshold:
                w1 = 1
            else:
                w1 = (1 - edgeTaper) * (1 - (math.fabs(length_in_direction1) / (max_length_in_dir1 / 2))) + edgeTaper
            if max_length_in_dir2 < threshold:
                w2 = 1
            else:
                w2 = (1 - edgeTaper) * (1 - (math.fabs(length_in_direction2) / (max_length_in_dir2 / 2))) + edgeTaper
        elif taper.lower() == 'hamming':  # Hamming Window
            if max_length_in_dir1 < threshold:
                w1 = 1
            else:
                w1 = 0.54 - 0.46 * math.cos(2 * math.pi * (length_in_direction1 / max_length_in_dir1 - 0.5))
            if max_length_in_dir2 < threshold:
                w2 = 1
            else:
                w2 = 0.54 - 0.46 * math.cos(2 * math.pi * (length_in_direction2 / max_length_in_dir2 - 0.5))
        else:
            return 0

        return w1 * w2

    def beamform(self, phi_scan=0, theta_scan=0):
        '''
        Returns far field pattern calculated for a specific phi/scan angle requested.
        This is calculated based on the lattice vector spacing and the embedded element
        patterns of a ca-ddm or fa-ddm array in HFSS.

        Parameters:


                phi_scan (float/int): spherical cs for desired scan angle of beam
                theta_scan (float/int): spherical cs for desired scan angle of beam

        Returns:
                all_qtys (dict): dictionary with reTheta, rePhi,RealizedGain, theta,phi output
        '''
        num_ports = len(self.all_port_names)
        self.FindArrayCenterAndEdge()

        c = 299792458
        k = (2 * math.pi * self.freq) / c

        # ---------------------- METHOD : CalculatePhaseShifts -------------------
        # Calculates phase shifts between array elements in A and B directions,
        # PhaseShiftA and PhaseShiftB, given Wave Vector (k), lattice vectors
        # (Ax, Ay, Bx, By), Scan angles (theta, phi) using formula below
        # Phase Shift A = - (Ax*k*sinθ*cosφ + Ay*k*sinθ*sinφ)
        # Phase Shift B = - (Bx*k*sinθ*cosφ + By*k*sinθ*sinφ)
        # ------------------------------------------------------------------------

        theta_scan = math.radians(theta_scan)
        phi_scan = math.radians(phi_scan)

        if self.lattice_vectors is not None:
            phase_shift_A_rad = -1 * ((self.Ax * k * math.sin(theta_scan) * math.cos(phi_scan))
                                      + (self.Ay * k * math.sin(theta_scan) * math.sin(phi_scan)))
            phase_shift_B_rad = -1 * ((self.Bx * k * math.sin(theta_scan) * math.cos(phi_scan))
                                      + (self.By * k * math.sin(theta_scan) * math.sin(phi_scan)))

        w_dict = {}
        w_dict_ang = {}
        w_dict_mag = {}
        array_positions = {}
        for port_name in self.all_port_names:
            index_str = get_array_index(port_name)
            a = int(index_str[0])
            b = int(index_str[1])

            w_mag = np.round(np.abs(self.AssignWeight(a, b, taper=self.taper)), 6)
            if self.lattice_vectors is None:
                pos = self.element_position_dict[port_name]
                w_ang = -1 * ((pos[0] * k * math.sin(theta_scan) * math.cos(phi_scan))
                              + (pos[1] * k * math.sin(theta_scan) * math.sin(phi_scan)))

            else:
                w_ang = (a * phase_shift_A_rad + b * phase_shift_B_rad)
            # ToDo check for driven modal or terminal
            w_dict[port_name] = np.sqrt(w_mag) * np.exp(1j * w_ang)
            w_dict_ang[port_name] = w_ang
            w_dict_mag[port_name] = w_mag
            if self.lattice_vectors is not None:
                array_positions[port_name] = self.ElementPosition(a, b)
            else:
                array_positions[port_name] = self.element_position_dict[port_name]

        length_of_ff_data = len(
            self.data_dict[self.all_port_names[0]]['rETheta'])  # check length of data by using first port

        rEtheta_fields = np.zeros((num_ports, length_of_ff_data), dtype=complex)
        rEphi_fields = np.zeros((num_ports, length_of_ff_data), dtype=complex)
        w = np.zeros((1, num_ports), dtype=complex)
        # create port mapping
        for n, port in enumerate(self.all_port_names):
            re_theta = self.data_dict[port]['rETheta']  # this is re_theta index of loaded data
            re_phi = self.data_dict[port]['rEPhi']  # this is re_ohi index of loaded data

            w[0][n] = w_dict[port]  # build 1xNumPorts array of weights

            rEtheta_fields[n] = re_theta
            rEphi_fields[n] = re_phi

            theta_range = self.data_dict[port]['Theta']
            phi_range = self.data_dict[port]['Phi']
            Ntheta = len(theta_range)
            Nphi = len(phi_range)

        rEtheta_fields_sum = np.dot(w, rEtheta_fields)
        rEtheta_fields_sum = np.reshape(rEtheta_fields_sum, (Ntheta, Nphi))

        rEphi_fields_sum = np.dot(w, rEphi_fields)
        rEphi_fields_sum = np.reshape(rEphi_fields_sum, (Ntheta, Nphi))

        self.all_qtys = {}
        self.all_qtys['rEPhi'] = rEphi_fields_sum
        self.all_qtys['rETheta'] = rEtheta_fields_sum
        self.all_qtys['rETotal'] = np.sqrt(
            np.power(np.abs(rEphi_fields_sum), 2) + np.power(np.abs(rEtheta_fields_sum), 2))
        self.all_qtys['Theta'] = theta_range
        self.all_qtys['Phi'] = phi_range
        self.all_qtys['nPhi'] = Nphi
        self.all_qtys['nTheta'] = Ntheta
        pin = np.sum(np.power(np.abs(w), 2))
        self.all_qtys['Pincident'] = pin
        print(f'Incident Power: {pin}')
        real_gain = 2 * np.pi * np.abs(np.power(self.all_qtys['rETotal'], 2)) / pin / 377
        self.all_qtys['RealizedGain'] = real_gain
        self.all_qtys['RealizedGain_dB'] = 10 * np.log10(real_gain)
        self.max_gain = np.max(10 * np.log10(real_gain))
        self.min_gain = np.min(10 * np.log10(real_gain))
        print(f'Peak Realized Gain: {self.max_gain} dB')
        self.all_qtys['Element_Location'] = array_positions

        return self.all_qtys

    def convert_units(value, oldUnits, newUnits):
        """
        used for converting between common unit types in HFSS
        """
        unitConv = {"nm": .000000001, "um": .000001, "mm": .001, "meter": 1.0, "cm": .01, "ft": .3048, "in": .0254, "mil": .0000254, "uin": .0000000254}
        value =float(value)
        sf = 1.0

        BaseUnits = None
        NewUnits = None
        if oldUnits.lower() in unitConv:
            BaseUnits = unitConv[oldUnits.lower()]
        if newUnits.lower() in unitConv:
            NewUnits = unitConv[newUnits.lower()]

        if BaseUnits != None and NewUnits != None:
            sf = BaseUnits/NewUnits


        if oldUnits != newUnits:

            nuValue = value*sf
        else:
            nuValue = value


        return nuValue