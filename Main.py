# -*- coding: utf-8 -*-
"""
Created on Tue Dec 14 14:48:07 2021

Script will generate the total far field patterns from multiple embedded element
patterns exported from HFSS. THe HFSS design needs to be setup using finite array
domain decomposition. The total far field pattern will be displayed, and the 
beam can interactivly be pointed in any direction, along with specific weighting
functions be defined. 

An eep file can also be exported for use in STK

@author: asligar
"""



from Lib.aedt_utils import AEDTutils
from Lib.FarField_Setup import FarField_Utils
from Lib.FarFieldsProcessing import Load_FF_Fields
from Lib.Reporter import Report_Module
from pyaedt import Hfss

###############################################################################
#
# BEGIN USER INPUTS
#
###############################################################################

export_eep = True # eep file used in STK
project_name = "5ghz_array"
design_name = "triangular_lattice"
solution_setup_name = 'Setup1 : LastAdaptive'
ff_setup_name = 'Infinite Sphere1'
freq = 5e9
taper='hamming' #can be cosine triangular hamming or flat

###############################################################################
#
# END USER INPUTS
#
###############################################################################

def main(project_name,design_name,solution_setup_name=None,ff_setup_name='Infinite_Sphere1',freq='',taper='flat',export_eep=False):


    with Hfss(non_graphical=False, new_desktop_session=False,specified_version='2022.1') as aedtapp:
        aedt = AEDTutils(aedtapp,project_name=project_name,design_name=design_name)

        ff_setup = FarField_Utils(aedt)
        if ff_setup_name not in aedt.ff_setups_names:
            ff_setup_name =ff_setup.insert_infinite_sphere(setup_name = ff_setup_name,overwrite = False,cs_name='Global')
        
        #get all teh embedded element
        eep_dict = ff_setup.export_all_ffd(ff_setup_name,
                                     freq=freq,
                                     setup_name = solution_setup_name,
                                     overwrite=False)
        
        lattice_vectors = ff_setup.get_lattice_vectors(aedt.model_units)
        
        scan_angle_theta = 0
        scan_angle_phi = 0
        ff = Load_FF_Fields(eep_dict,lattice_vectors)
        ff.taper = taper
        ff.freq = freq
        results = ff.beamform(phi_scan=scan_angle_phi,theta_scan=scan_angle_theta)
        
        reports = Report_Module(aedt,results,ff.data_dict)
        # reports.plot_far_field_rect(qty_str='RealizedGain',convert_to_dB=True)
        
        # reports.plot_2d_cut(phi='all',theta=scan_angle_theta,
        #                 qty_str = 'RealizedGain',
        #                 title='Azimuth', 
        #                 convert_to_dB=True)

        # reports.plot_2d_cut(phi=scan_angle_phi,theta='all',
        #                 qty_str = 'RealizedGain',
        #                 title='Elevation', 
        #                 convert_to_dB=True)
        
        # reports.polar_plot_3d(qty_str = 'RealizedGain',
        #                     convert_to_dB=True)
        
        cad_mesh = reports.get_geometry()
        reports.polar_plot_3d_pyvista(ff,qty_str = 'RealizedGain',
                                            convert_to_dB=True,
                                            cad_mesh=cad_mesh)
        
        
        #I added a simple demo to demonstrate 2 beams being excited at the same time
        #uncomment the lines below to try the demo
        
        # reports.polar_plot_3d_pyvista_2beams(ff,qty_str = 'RealizedGain',
        #                                     convert_to_dB=True,
        #                                     cad_mesh=cad_mesh)
        
        if export_eep:
            path_to_eep = ff.save_eep(aedt.aedtapp.project_path,aedt.aedtapp.design_name)


if __name__ == "__main__":
    main(project_name,
         design_name,
         solution_setup_name=solution_setup_name,
         ff_setup_name=ff_setup_name,
         freq=freq,
         taper=taper,
         export_eep = export_eep)