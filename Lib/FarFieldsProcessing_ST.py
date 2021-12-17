# -*- coding: utf-8 -*-
"""
Created on Tue Sep  7 15:50:47 2021

@author: asligar
"""
import numpy as np
import time as walltime
import os
import math




def load(ffd_dict,lattice_vectors=[-1, 1, 0,-1,-1, 0.0]):
    
    lattice_vectors=lattice_vectors
    
    Ax = float(lattice_vectors[0])
    Ay = float(lattice_vectors[1])
    Bx = float(lattice_vectors[3])
    By = float(lattice_vectors[4])
    
    ffd_dict = ffd_dict
    data_dict = {}
    all_port_names = list(ffd_dict.keys())
    
    valid_ffd = True
    all_ports = list(ffd_dict.keys())
    if os.path.exists(ffd_dict[all_ports[0]]):
        with open(ffd_dict[all_ports[0]], 'r') as reader:
            theta=[int(i) for i in reader.readline().split()] 
            phi=[int(i) for i in reader.readline().split()]
            num_freq=int(reader.readline().split()[1])
            frequency=float(reader.readline().split()[1])
        reader.close()
        results_dict = {}
        for port in ffd_dict.keys():
            if ':' in port:
                port = port.split(':')[0]
            temp_dict = {}
            theta_range=np.linspace(*theta)
            phi_range= np.linspace(*phi)
            
            ntheta=len(theta_range)
            nphi=len(phi_range)
            
            if os.path.exists(ffd_dict[port]):
                eep_txt=np.loadtxt(ffd_dict[port], skiprows=4)
                Etheta=np.vectorize(complex)(eep_txt[:,0], eep_txt[:,1])
                Ephi=np.vectorize(complex)(eep_txt[:,2], eep_txt[:,3])
                
                #eep=np.column_stack((etheta, ephi))  
                
                temp_dict['Theta']=theta_range
                temp_dict['Phi'] = phi_range
                temp_dict['rETheta']=Etheta
                temp_dict['rEPhi']=Ephi
                
                data_dict[port]=temp_dict
            else:
                valid_ffd=False
            
        #differential area of sphere, based on observation angle
        d_theta = np.abs(theta_range[1]-theta_range[0])
        d_phi = np.abs(phi_range[1]-phi_range[0])
        diff_area=np.radians(d_theta)*np.radians(d_phi)*np.sin(np.radians(theta_range)) 
        num_samples = len(temp_dict['rETheta'])
        
        return data_dict
    else:
        valid_ffd=False
        print('ERROR: Far Field Files are Missing')
        return False





def GetArrayIndex(port_name):
        str1 = port_name.split('[', 1)[1].split(']', 1)[0]
        index_str = str1.split(',')
        return index_str
def FindArrayMinMax( all_port_names=None):
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
        all_port_names = all_port_names
    row_min = 1
    row_max = 1
    col_min = 1
    col_max = 1
    rows = []
    cols = []
    for portstring in all_port_names:
        index_str = GetArrayIndex(portstring)
        rows.append( int(index_str[1]) )
        cols.append(int(index_str[0]))

    row_min = np.min(rows)
    row_max = np.max(rows)
    col_min = np.min(cols)
    col_max = np.max(cols)
    return [col_min, col_max, row_min,row_max]

def FindArrayCenterAndEdge(all_port_names,Ax,Ay,Bx,By):
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


    for i in range(0, len(all_port_names)):                           
        index_str = GetArrayIndex(all_port_names[i])
        row = int(index_str[1]) 
        col = int(index_str[0])
        a = row
        b = col

        activeCells.append((a,b)) #because ffd is assuming all ffd files are active
    if len(activeCells) == 0: return

    [a_min,a_max,b_min,b_max] = FindArrayMinMax(all_port_names)
    
    CenterA = (a_min+a_max)/2
    CenterB = (b_min+b_max)/2
    CenterX = (CenterA+0.5) * Ax + (CenterB+0.5) * Bx
    CenterY = (CenterA+0.5) * Ay + (CenterB+0.5) * By

    CenterA = CenterA
    CenterB = CenterB
    CenterX = CenterX
    CenterY = CenterY
    # find the distance from the edge to the center
    AMax = a_max-a_min
    BMax = b_max-b_min
    
    AMax = AMax
    BMax = BMax
    for a,b in activeCells:
        x = (a + 0.5) * Ax + (b + 0.5) * Bx
        y = (a + 0.5) * Ay + (b + 0.5) * By
        x_dis = abs(x-CenterX)
        y_dis = abs(y-CenterY)
        distance = math.sqrt(x_dis**2 + y_dis**2)
        XMax = max(XMax,x_dis)
        YMax = max(YMax,y_dis)
        RMax = max(RMax,distance)

    RMax = RMax
    XMax = XMax
    YMax = YMax
    RMax *= 2
    XMax *= 2
    YMax *= 2

    return CenterA, CenterB, CenterX, CenterY,AMax,BMax
def ElementPosition( a, b,Ax,Ay,Bx,By,CenterX,CenterY):
    a = int(a)
    b = int(b)


    x = (a + 0.5) * Ax + (b + 0.5) * Bx
    y = (a + 0.5) * Ay + (b + 0.5) * By        
    x_dis = x-CenterX
    y_dis = y-CenterY
    
    return np.array([x_dis,y_dis,0])
def AssignWeight( a, b,CenterA,CenterB,AMax,BMax,taper='flat'):
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
    length_in_direction1 = a - CenterA
    length_in_direction2 = b - CenterB
    max_length_in_dir1 = AMax
    max_length_in_dir2 = BMax 

    
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

def beamform(data_dict,lattice_vectors,freq,phi_scan=0,theta_scan=0,taper='cosine'):
    '''
    Returns far field pattern calculated for a specific phi/scan angle requested.
    This is calculated based on the lattice vector spacing and the embedded element
    patterns of a ca-ddm or fa-ddm array in HFSS.

    Parameters:
            ff_data (dict): dictionary of embedded element patterns
            freq (float): frequency to calculate far field for
            phi_scan (float/int): spherical cs for desired scan angle of beam
            theta_scan (float/int): spherical cs for desired scan angle of beam
            taper [str]: type of aperture taper to be applied
    Returns:
            all_qtys (dict): dictionary with reTheta, rePhi,RealizedGain, theta,phi output
    '''
    num_ports = len(data_dict)
    all_port_names = list(data_dict.keys())
    
    Ax = float(lattice_vectors[0])
    Ay = float(lattice_vectors[1])
    Bx = float(lattice_vectors[3])
    By = float(lattice_vectors[4])
    
    CenterA, CenterB, CenterX, CenterY, AMax, BMax = FindArrayCenterAndEdge(all_port_names,Ax,Ay,Bx,By)

    c=299792458
    k = (2*math.pi*freq)/c

    #---------------------- METHOD : CalculatePhaseShifts -------------------
    # Calculates phase shifts between array elements in A and B directions,
    # PhaseShiftA and PhaseShiftB, given Wave Vector (k), lattice vectors
    # (Ax, Ay, Bx, By), Scan angles (theta, phi) using formula below
    # Phase Shift A = - (Ax*k*sinθ*cosφ + Ay*k*sinθ*sinφ)
    # Phase Shift B = - (Bx*k*sinθ*cosφ + By*k*sinθ*sinφ)
    #------------------------------------------------------------------------
    
    
    theta_scan = math.radians(theta_scan)
    phi_scan = math.radians(phi_scan)

    phase_shift_A_rad = -1*( (Ax*k*math.sin(theta_scan)*math.cos(phi_scan)) 
                             + (Ay*k*math.sin(theta_scan)*math.sin(phi_scan)) )
    phase_shift_B_rad = -1*( (Bx*k*math.sin(theta_scan)*math.cos(phi_scan)) 
                             + (By*k*math.sin(theta_scan)*math.sin(phi_scan)) )
    
   
    w_dict ={}
    w_dict_ang = {}
    w_dict_mag = {}
    array_positions = {}
    for port_name in all_port_names:

        index_str = GetArrayIndex(port_name)
        a = int(index_str[0])
        b = int(index_str[1])
        w_mag = np.round(np.abs(AssignWeight(a, b,CenterA,CenterB,AMax, BMax,taper=taper)),3)
        w_ang = (a*phase_shift_A_rad+b*phase_shift_B_rad)
        # ToDo check for driven modal or terminal
        w_dict[port_name] = np.sqrt(w_mag)*np.exp(1j*w_ang)
        w_dict_ang[port_name] = w_ang
        w_dict_mag[port_name] = w_mag
        array_positions[port_name] = ElementPosition(a,b,Ax,Ay,Bx,By,CenterX, CenterY)

    #print(list(data_dict.keys()))
    length_of_ff_data = len(data_dict[all_port_names[0]]['rETheta']) #check lenght of data by usingfirst port
    
    rEtheta_fields= np.zeros((num_ports,length_of_ff_data),dtype=complex)
    rEphi_fields= np.zeros((num_ports,length_of_ff_data),dtype=complex)
    w= np.zeros((1,num_ports),dtype=complex)
    #create port mapping
    for n, port in enumerate(all_port_names):
        re_theta = data_dict[port]['rETheta'] #this is re_theta index of loaded data
        re_phi = data_dict[port]['rEPhi'] #this is re_ohi index of loaded data
    
        w[0][n]=w_dict[port] #build 1xNumPorts array of weights
    
        rEtheta_fields[n] = re_theta
        rEphi_fields[n] = re_phi
    
        theta_range=data_dict[port]['Theta']
        phi_range=data_dict[port]['Phi']
        Ntheta=len(theta_range)
        Nphi=len(phi_range)
    
    rEtheta_fields_sum = np.dot(w,rEtheta_fields)
    rEtheta_fields_sum  = np.reshape(rEtheta_fields_sum ,(Ntheta,Nphi))

    rEphi_fields_sum  = np.dot(w,rEphi_fields)
    rEphi_fields_sum  = np.reshape(rEphi_fields_sum ,(Ntheta,Nphi))
    
    all_qtys={}
    #all_qtys['rEPhi'] = rEphi_fields_sum 
    #all_qtys['rETheta'] = rEtheta_fields_sum 
    all_qtys['rETotal'] = np.sqrt(np.power(np.abs(rEphi_fields_sum ),2)+np.power(np.abs(rEtheta_fields_sum ),2))
    all_qtys['Theta'] = theta_range
    all_qtys['Phi'] = phi_range
    all_qtys['nPhi'] = Nphi
    all_qtys['nTheta'] = Ntheta
    pin=np.sum(np.power(np.abs(w),2))
    all_qtys['Pincident'] = pin
    print(f'Incident Power: {pin}')
    real_gain = 2*np.pi*np.abs(np.power(all_qtys['rETotal'],2))/pin/377
    all_qtys['RealizedGain'] = real_gain
    #all_qtys['RealizedGain_dB'] = 10*np.log10(real_gain)
    all_qtys['Element_Location'] = array_positions
    
    return all_qtys