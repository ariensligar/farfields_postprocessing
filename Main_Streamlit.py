# -*- coding: utf-8 -*-
"""
Created on Tue Dec 14 14:48:07 2021

Script will generate terrain and import shape file into HFSS

@author: asligar
"""

from matplotlib.backends.backend_agg import RendererAgg
import streamlit as st

from Lib.aedt_utils import AEDTutils
from Lib.FarField_Setup import FarField_Utils
import Lib.FarFieldsProcessing_ST as ff
import Lib.Geo_ST as geo
from Lib.Reporter_ST import Report_Module
from pyaedt import Hfss
import streamlit.components.v1 as components
from itkwidgets import view
from ipywidgets import embed
import pyvista as pv
import vtk
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

###############################################################################
#
# END USER INPUTS
#
###############################################################################

#st.set_page_config(layout="wide")
st.title("Far Field Post Processing - AEDT")
st.markdown(
"""
This is a test to create interactive plots in AEDT using Streamlit
""")
global aedtapp

@st.cache()
def get_eep(ff_setup_name):
    def get_element_positions(temp_dict,lattice_vectors):
        array_positions= {}
        all_port_names = list(temp_dict.keys())
        
        Ax = float(lattice_vectors[0])
        Ay = float(lattice_vectors[1])
        Bx = float(lattice_vectors[3])
        By = float(lattice_vectors[4])
        
        CenterA, CenterB, CenterX, CenterY, AMax, BMax = ff.FindArrayCenterAndEdge(all_port_names,Ax,Ay,Bx,By)
        
        for port_name in all_port_names:
            index_str = ff.GetArrayIndex(port_name)
            a = int(index_str[0])
            b = int(index_str[1])
            array_positions[port_name] = ff.ElementPosition(a,b,Ax,Ay,Bx,By,CenterX, CenterY)
        return array_positions
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
        array_positions = get_element_positions(eep_dict, lattice_vectors)
        meshes, all_max = geo.get_geometry(aedt,array_positions)
        
        return eep_dict, lattice_vectors, meshes, all_max




@st.cache()
def load_fields(temp_dict, lattice_vectors):
    data_dict = ff.load(temp_dict,lattice_vectors)
    
    return data_dict



phi=st.sidebar.slider('Phi',min_value=-180,max_value=180,value=0)
theta=st.sidebar.slider('Theta',min_value=-180,max_value=180,value=0)
#st.write(x,'asdfasdf',min)

eep_dict, lattice_vectors, meshes, all_max = get_eep(ff_setup_name)


data  = load_fields(eep_dict, lattice_vectors )




results = ff.beamform(data,lattice_vectors,freq,phi_scan=phi,theta_scan=theta,taper='flat')



reports = Report_Module(results)
# reports.plot_far_field_rect(qty_str='RealizedGain',convert_to_dB=True)
fig, ax = reports.polar_plot_3d(qty_str = 'RealizedGain',
                    convert_to_dB=True)
st.pyplot(fig)
fig, ax = reports.plot_2d_cut(phi='all',theta=theta,
                qty_str = 'RealizedGain',
                title='Azimuth', 
                convert_to_dB=True)
st.pyplot(fig)

array_cad = pv.read(meshes)




#     # fig, ax = reports.plot_2d_cut(phi=phi,theta='all',
#     #                 qty_str = 'RealizedGain',
#     #                 title='Elevation', 
#     #                 convert_to_dB=True)
#     # st.pyplot(fig)


cad_and_ff,mesh = reports.polar_plot_3d_pyvista(qty_str = 'RealizedGain',
                                    convert_to_dB=True,
                                    cad_mesh=array_cad)

    
view_width = 800
view_height = 600

snippet = embed.embed_snippet(views=view(mesh))
html = embed.html_template.format(title="", snippet=snippet)
components.html(html, width=view_width, height=view_height)