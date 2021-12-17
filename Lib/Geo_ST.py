# -*- coding: utf-8 -*-
"""
Created on Thu Dec 16 08:53:09 2021

@author: asligar
"""

import os
import Lib.Utillities as utils
import pyvista as pv
import numpy as np

def get_geometry(aedt,array_positions,is_antenna_array=True):
    oEditor = aedt.aedtapp.odesign.SetActiveEditor("3D Modeler")
    
    #obj is being exported as model units, scaling factor needed for display
    model_units = oEditor.GetModelUnits()
    sf = utils.convert_units(1,'meter',model_units)

    
    bounding_box = oEditor.GetModelBoundingBox()
    xmax = float(bounding_box[3])-float(bounding_box[0])
    ymax = float(bounding_box[4])-float(bounding_box[1])
    zmax = float(bounding_box[5])-float(bounding_box[2])

    
    geo_path = aedt.aedtapp.results_directory+'/geo/'
    if not os.path.exists(geo_path):
        os.makedirs(geo_path)
    cad_file = geo_path + '/geometry.obj'
    
    #get list of object we want to display
    non_model_objects = oEditor.GetObjectsInGroup('Non Model')
    all_objects = oEditor.GetMatchedObjectName('*')
    
    s = set(non_model_objects)
    model_objects = [x for x in all_objects if x not in s]
    air_objects = oEditor.GetObjectsByMaterial('vacuum')
    air_objects += oEditor.GetObjectsByMaterial('air')
    
    #don't display objects that are vacuum or air
    reduced_model_objects = model_objects
    for each in air_objects:
        if each in model_objects: 
            reduced_model_objects.remove(each)
            
    objects_to_display = []
    selected_objects = oEditor.GetSelections()
    print('TIP: Geometry selected in AEDT will be displayed along with far field pattern')
    print('TIP: IF no selected geometry, all model objects will be displayed')
    if len(selected_objects)>=1:
        objects_to_display = selected_objects
    else:
        objects_to_display = reduced_model_objects
    print("INFO: Exporting Geometry for Display")
    oEditor.ExportModelMeshToFile(cad_file, objects_to_display)
    print("...Done")
    
    meshes= pv.PolyData()
    if os.path.exists(cad_file):
        cad_mesh = pv.read(cad_file)
        color_display_type = ''
        if is_antenna_array:
            for each in array_positions:
                translated_mesh=cad_mesh.copy()
                offset_xyz = array_positions[each]*sf
                if np.abs(2*offset_xyz[0])>xmax:#assume array is centere, factor of 2
                    xmax=offset_xyz[0]*2
                if np.abs(2*offset_xyz[1])>ymax:#assume array is centere, factor of 2
                    ymax=offset_xyz[1]*2
                translated_mesh.translate(offset_xyz)
                meshes+=translated_mesh
        else:
                translated_mesh.translate(array_positions[each]*sf)
                meshes=translated_mesh
    else:
        print('WARNING: Unable to display CAD Geometry, ' + cad_file + ' is not found')
    
    output_geo = geo_path + '/geometry_full.vtk'
    meshes.save(output_geo,binary=True)
    all_max = np.max(np.array([xmax,ymax,zmax]))
    return output_geo, all_max