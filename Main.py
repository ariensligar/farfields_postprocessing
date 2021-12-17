# -*- coding: utf-8 -*-
"""
Created on Tue Dec 14 14:48:07 2021

Script will generate terrain and import shape file into HFSS

@author: asligar
"""



from Lib.aedt_utils import AEDTutils
from Lib.FarField_Setup import FarField_Utils
from Lib.FarFieldsProcessing import Load_FF_Fields
from Lib.Reporter import Report_Module
from pyaedt import Hfss

import pyvista as pv

###############################################################################
#
# BEGIN USER INPUTS
#
###############################################################################


project_name = "5ghz_array"
design_name = "triangular_lattice_small"
solution_setup_name = 'Setup1 : LastAdaptive'
ff_setup_name = 'Infinite Sphere1'
freq = 5e9
taper='hamming' #can be cosine triangular hamming or flat

###############################################################################
#
# END USER INPUTS
#
###############################################################################

def main(project_name,design_name,solution_setup_name=None,ff_setup_name='Infinite_Sphere1',freq='',taper='flat'):


    with Hfss(non_graphical=False, new_desktop_session=False,specified_version='2022.1') as aedtapp:
        aedt = AEDTutils(aedtapp,project_name=project_name,design_name=design_name)

        ff_setup = FarField_Utils(aedt)
        if ff_setup_name not in aedt.ff_setups_names:
            ff_setup_name =ff_setup.insert_infinite_sphere(setup_name = ff_setup_name,overwrite = False,cs_name='Global')
        
        eep_dict = ff_setup.export_all_ffd(ff_setup_name,
                                     freq=freq,
                                     setup_name = solution_setup_name,
                                     overwrite=False)
        
        lattice_vectors = ff_setup.get_lattice_vectors(aedt.model_units)
        
        scan_angle_theta = 0
        scan_angle_phi = 0
        ff = Load_FF_Fields(eep_dict,lattice_vectors)
        ff.taper='flat'
        ff.freq = freq
        results = ff.beamform(phi_scan=scan_angle_phi,theta_scan=scan_angle_theta)
        
        reports = Report_Module(aedt,results,ff.data_dict)
        reports.plot_far_field_rect(qty_str='RealizedGain',convert_to_dB=True)
        
        reports.plot_2d_cut(phi='all',theta=scan_angle_theta,
                        qty_str = 'RealizedGain',
                        title='Azimuth', 
                        convert_to_dB=True)

        reports.plot_2d_cut(phi=scan_angle_phi,theta='all',
                        qty_str = 'RealizedGain',
                        title='Elevation', 
                        convert_to_dB=True)
        
        reports.polar_plot_3d(qty_str = 'RealizedGain',
                            convert_to_dB=True)
        
        cad_mesh = reports.get_geometry()
        cad_and_ff = reports.polar_plot_3d_pyvista(ff,qty_str = 'RealizedGain',
                                            convert_to_dB=True,
                                            cad_mesh=cad_mesh)
        #cad_and_ff.show()
        
if __name__ == "__main__":
    main(project_name,
         design_name,
         solution_setup_name=solution_setup_name,
         ff_setup_name=ff_setup_name,
         freq=freq,
         taper=taper)