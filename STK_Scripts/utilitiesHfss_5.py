# -*- coding: utf-8 -*-
"""

Contains function for extracting far field values and embedded element patterns

Arien Sligar (arien.sligar@ansys.com)
Last Update: 06/24/2020
"""
import numpy as np
import os
import math
import time as walltime
import pickle


class utils():
    def __init__(self,all_eep,lattice_vectors=[-1, 1, 0,-1,-1, 0.0]):
        '''
        Loads all the near fields into memory.
        Parameters
        ----------
        nfd_files_dict : dict
            dictionary of all near field files, with key=port name.
        length : float
            size of near field rectangle, this will coorponds to the x-axis. of
            the local CS used in HFSS to reference the near field surface
        width : float
            size of near field rectangle, this will coorponds to the y-axis. of
            the local CS used in HFSS to reference the near field surface
        length_n : int
            number of points in the near field rectangle.
        width_n : int
            number of points in the near field rectangle..
        Returns
        -------
        None.
        creates a dictionary of (data_dict) of near field values and positions
        for each port
        '''
        time_before = walltime.time()
        print('Loading Embedded Element Patterns...')
        self.data_for_eep_export = {}
        self.freq = 1e9
        self.taper='flat'
        valid_ffd = True
        all_ports = list(all_eep.keys())
        #differential area of sphere, based on observation angle
        self.all_port_names = list(all_eep.keys())
        self.solution_type = 'DrivenModal'
        self.unique_beams = None
        self.renormalize = False
        self.renormalize_dB = True
        self.renorm_value= 1
        self.data_dict = {}
        for port in all_eep.keys():
            if ':' in port:
                port = port.split(':')[0]
            temp_dict = {}
            temp_dict['Theta']=all_eep[port][0]
            temp_dict['Phi'] = all_eep[port][1]
            temp_dict['rETheta']=all_eep[port][2]
            temp_dict['rEPhi']=all_eep[port][3]
            self.data_dict[port]=temp_dict
        self.valid_ffd = valid_ffd
        self.lattice_vectors = lattice_vectors
        self.Ax = float(lattice_vectors[0])
        self.Ay = float(lattice_vectors[1])
        self.Bx = float(lattice_vectors[3])
        self.By = float(lattice_vectors[4])
        elapsed_time = walltime.time()-time_before
        print(f'Loading Embedded Element Patterns...Done: {elapsed_time}seconds')
    def GetArrayIndex(self,port_name):
            str1 = port_name.split('[', 1)[1].split(']', 1)[0]
            index_str = str1.split(',')
            return index_str
    def FindArrayMinMax(self, all_port_names=None):
        '''
        Parameters
        ----------
        all_port_names : list of strings, if this isn't suplied we will assume
        it was already intialized in the main script before calling ff_beamsteer()
            DESCRIPTION. The default is None.
        Returns
        -------
        [int, int, int, int]
            [column min, column max, row min, row max] of array
        '''
        if all_port_names==None:
            all_port_names = self.all_port_names
        row_min = 1
        row_max = 1
        col_min = 1
        col_max = 1
        rows = []
        cols = []
        for portstring in all_port_names:
            index_str = self.GetArrayIndex(portstring)
            rows.append( int(index_str[1]) )
            cols.append(int(index_str[0]))
        row_min = np.min(rows)
        row_max = np.max(rows)
        col_min = np.min(cols)
        col_max = np.max(cols)
        return [col_min, col_max, row_min,row_max]
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
            index_str = self.GetArrayIndex(self.all_port_names[i])
            row = int(index_str[1]) 
            col = int(index_str[0])
            a = row
            b = col
            activeCells.append((a,b)) #because ffd is assuming all ffd files are active
        if len(activeCells) == 0: return
        [a_min,a_max,b_min,b_max] = self.FindArrayMinMax()
        CenterA = (a_min+a_max)/2
        CenterB = (b_min+b_max)/2
        CenterX = (CenterA+0.5) * self.Ax + (CenterB+0.5) * self.Bx
        CenterY = (CenterA+0.5) * self.Ay + (CenterB+0.5) * self.By
        self.CenterA = CenterA
        self.CenterB = CenterB
        self.CenterX = CenterX
        self.CenterY = CenterY
        # find the distance from the edge to the center
        AMax = a_max-a_min
        BMax = b_max-b_min
        self.AMax = AMax
        self.BMax = BMax
        for a,b in activeCells:
            x = (a + 0.5) * self.Ax + (b + 0.5) * self.Bx
            y = (a + 0.5) * self.Ay + (b + 0.5) * self.By
            x_dis = abs(x-CenterX)
            y_dis = abs(y-CenterY)
            distance = math.sqrt(x_dis**2 + y_dis**2)
            XMax = max(XMax,x_dis)
            YMax = max(YMax,y_dis)
            RMax = max(RMax,distance)
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
        x_dis = x-self.CenterX
        y_dis = y-self.CenterY
        return np.array([x_dis,y_dis,0])
    def AssignWeight(self, a, b,taper='flat'):
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
        if taper.lower() == 'flat':                                  #Flat
            return 1
        cosinePow = 1
        edgeTaper_dB=-200
        edgeTaper = 10**((float(edgeTaper_dB))/20)
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
        if taper.lower() == 'cosine':                                   #Cosine
            if max_length_in_dir1 < threshold: w1 = 1
            else: w1 = (1-edgeTaper)*(math.cos(math.pi*length_in_direction1/max_length_in_dir1))**cosinePow + edgeTaper
            if max_length_in_dir2 < threshold: w2 = 1
            else: w2 = (1-edgeTaper)*(math.cos(math.pi*length_in_direction2/max_length_in_dir2))**cosinePow + edgeTaper
        elif taper.lower() == 'triangular':                                #Triangular
            if max_length_in_dir1 < threshold: w1 = 1
            else: w1 = (1-edgeTaper)*(1-(math.fabs(length_in_direction1)/(max_length_in_dir1/2))) + edgeTaper
            if max_length_in_dir2 < threshold: w2 = 1
            else: w2 = (1-edgeTaper)*(1-(math.fabs(length_in_direction2)/(max_length_in_dir2/2))) + edgeTaper
        elif taper.lower() == 'hamming':                                #Hamming Window
                if max_length_in_dir1 < threshold: w1 = 1
                else: w1 = 0.54 - 0.46 * math.cos(2*math.pi*(length_in_direction1/max_length_in_dir1-0.5))
                if max_length_in_dir2 < threshold: w2 = 1
                else: w2 = 0.54 - 0.46 * math.cos(2*math.pi*(length_in_direction2/max_length_in_dir2-0.5))
        else:
            return 0
        return w1*w2
    def beamform(self,phi_scan=0,theta_scan=0):
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
        c=299792458
        k = (2*math.pi*self.freq)/c
        #---------------------- METHOD : CalculatePhaseShifts -------------------
        # Calculates phase shifts between array elements in A and B directions,
        # PhaseShiftA and PhaseShiftB, given Wave Vector (k), lattice vectors
        # (Ax, Ay, Bx, By), Scan angles (theta, phi) using formula below
        # Phase Shift A = - (Ax*k*sinθ*cosφ + Ay*k*sinθ*sinφ)
        # Phase Shift B = - (Bx*k*sinθ*cosφ + By*k*sinθ*sinφ)
        #------------------------------------------------------------------------
        theta_scan = math.radians(theta_scan)
        phi_scan = math.radians(phi_scan)
        phase_shift_A_rad = -1*( (self.Ax*k*math.sin(theta_scan)*math.cos(phi_scan)) 
                                 + (self.Ay*k*math.sin(theta_scan)*math.sin(phi_scan)) )
        phase_shift_B_rad = -1*( (self.Bx*k*math.sin(theta_scan)*math.cos(phi_scan)) 
                                 + (self.By*k*math.sin(theta_scan)*math.sin(phi_scan)) )
        w_dict ={}
        w_dict_ang = {}
        w_dict_mag = {}
        array_positions = {}
        for port_name in self.all_port_names:
            index_str = self.GetArrayIndex(port_name)
            a = int(index_str[0])
            b = int(index_str[1])
            w_mag = np.round(np.abs(self.AssignWeight(a, b,taper=self.taper)),6)
            w_ang = (a*phase_shift_A_rad+b*phase_shift_B_rad)
            # ToDo check for driven modal or terminal
            w_dict[port_name] = np.sqrt(w_mag)*np.exp(1j*w_ang)
            w_dict_ang[port_name] = w_ang
            w_dict_mag[port_name] = w_mag
            array_positions[port_name] = self.ElementPosition(a,b)
        length_of_ff_data = len(self.data_dict[self.all_port_names[0]]['rETheta']) #check lenght of data by usingfirst port
        rEtheta_fields= np.zeros((num_ports,length_of_ff_data),dtype=complex)
        rEphi_fields= np.zeros((num_ports,length_of_ff_data),dtype=complex)
        w= np.zeros((1,num_ports),dtype=complex)
        #create port mapping
        for n, port in enumerate(self.all_port_names):
            re_theta = self.data_dict[port]['rETheta'] #this is re_theta index of loaded data
            re_phi = self.data_dict[port]['rEPhi'] #this is re_ohi index of loaded data
            w[0][n]=w_dict[port] #build 1xNumPorts array of weights
            rEtheta_fields[n] = re_theta
            rEphi_fields[n] = re_phi
            theta_range=self.data_dict[port]['Theta']
            phi_range=self.data_dict[port]['Phi']
            Ntheta=len(theta_range)
            Nphi=len(phi_range)
        rEtheta_fields_sum = np.dot(w,rEtheta_fields)
        rEtheta_fields_sum  = np.reshape(rEtheta_fields_sum ,(Ntheta, Nphi))
        rEphi_fields_sum  = np.dot(w,rEphi_fields)
        rEphi_fields_sum  = np.reshape(rEphi_fields_sum ,(Ntheta, Nphi))
        self.all_qtys={}
        self.all_qtys['rEPhi'] = rEphi_fields_sum
        self.all_qtys['rETheta'] = rEtheta_fields_sum
        self.all_qtys['rETotal'] = np.sqrt(np.power(np.abs(rEphi_fields_sum ),2)+np.power(np.abs(rEtheta_fields_sum ),2))
        self.all_qtys['Theta'] = theta_range
        self.all_qtys['Phi'] = phi_range
        self.all_qtys['nPhi'] = Nphi
        self.all_qtys['nTheta'] = Ntheta
        pin=np.sum(np.power(np.abs(w),2))
        self.all_qtys['Pincident'] = pin
        print(f'Incident Power: {pin}')
        real_gain = 2*np.pi*np.abs(np.power(self.all_qtys['rETotal'],2))/pin/377
        self.all_qtys['RealizedGain'] = real_gain
        self.all_qtys['RealizedGain_dB'] = 10*np.log10(real_gain)
        self.max_gain = np.max(10*np.log10(real_gain))
        self.min_gain = np.min(10*np.log10(real_gain))
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