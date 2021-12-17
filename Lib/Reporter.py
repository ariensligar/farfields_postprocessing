# -*- coding: utf-8 -*-
"""
Created on Tue Sep  7 15:50:47 2021

@author: asligar
"""
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors

from matplotlib import cm
import pyvista as pv
import math

import Lib.Utillities as utils
import time as walltime

class Report_Module():
    def __init__(self,aedt,data,eep_data):
        self.aedtapp= aedt.aedtapp
        self.levels = 64
        self.data = data
        self.all_max =1

    def max_vs_beam_line(self,pd_max,title='Max Power Density',
                            pd_type_label = 'PD', 
                            save_name ="max_pd_line",
                            save_plot=False,
                            show_plot=True):
        beam_ids = list(pd_max.keys())
        pd_max_vals = list(pd_max.values())


        fig, ax = plt.subplots()
        ax.plot(beam_ids,pd_max_vals)
        
        ax.set(xlabel='Beam IDs', ylabel=pd_type_label,
               title=title)
        ax.grid()
        
        if save_plot:
            save_name_full = self.full_path + save_name + '.png'
            save_name_relative = self.relative_path + save_name + '.png'
            plt.savefig(save_name_full,dpi=300)
            self.all_figure_paths.append(save_name_relative)
        if show_plot:
            plt.show()
            
        if not show_plot:
            plt.close('all')

 
            
    def plot_far_field_rect(self,qty_str='RealizedGain',title='RectangularPlot',convert_to_dB=True):
        
        fig, ax = plt.subplots(figsize=(5, 5))
        
        if qty_str=='':
            qty_to_plot = self.data
            qty_str = 'Data'
        else:
            qty_to_plot = self.data[qty_str]
        qty_to_plot = np.reshape(qty_to_plot,(self.data['nTheta'],self.data['nPhi']))
        th,ph = np.meshgrid(self.data['Theta'], self.data['Phi'])

        if convert_to_dB:
            factor =20
            if 'Gain' in qty_str:
                factor =10
            qty_to_plot = factor*np.log10(np.abs(qty_to_plot))

        if title=='':            
            plt.title(qty_str)
        else:
            plt.title(title)
            
        plt.xlabel('Theta (degree)')
        plt.ylabel('Phi (degree)')

        plt.contourf(th,ph,qty_to_plot.T,levels=self.levels ,cmap='jet',)
        
        plt.colorbar()

        print('Peak '+ qty_str + ': ' +  str(np.max(qty_to_plot)))


    def plot_2d_cut(self,phi=0,theta='all',
                    qty_str = 'RealizedGain',
                    title='Far Field Cut', 
                    convert_to_dB=True):


        data_to_plot = self.data[qty_str]
        
        
        if phi=='all':
            
            x = self.data['Phi']
            theta_idx = self.find_nearest(self.data['Theta'],theta)
            y=data_to_plot[theta_idx]
            xlabel = 'Phi'
        if theta=='all':
            phi_idx = self.find_nearest(self.data['Phi'],phi)
            x = self.data['Theta']
            temp =data_to_plot.T
            y=temp[phi_idx]
            xlabel = 'Theta'
        
        
        if convert_to_dB:
            y=10*np.log10(y)
        fig, ax = plt.subplots()
        ax.plot(x,y)
    
        ax.set(xlabel=xlabel, ylabel=qty_str,
               title=title)
        ax.grid()
        

        plt.show()
        return fig, ax
    
    def plot_xy(self,x,y,title='xy plot', 
                            xlabel = 'x',
                            ylabel= 'y',
                            convert_to_dB=True):


        if convert_to_dB:
            x=10*np.log10(x)
        fig, ax = plt.subplots()
        ax.plot(x,y)
    
        ax.set(xlabel=xlabel, ylabel=ylabel,
               title=title)
        ax.grid()
        

        plt.show()


    def polar_plot_3d(self,qty_str = 'RealizedGain',
                        convert_to_dB=True):
    
        if convert_to_dB:
            ff_data = 10*np.log10(self.data[qty_str])
            #renormalize to 0 and 1 
            ff_max_dB = np.max(ff_data)
            ff_min_dB = np.min(ff_data)
            ff_data_renorm = (ff_data-ff_min_dB)/(ff_max_dB-ff_min_dB)
        else:
            ff_data = self.data[qty_str]
            #renormalize to 0 and 1 
            ff_max = np.max(ff_data)
            ff_min = np.min(ff_data)
            ff_data_renorm = (ff_data-ff_max)/(ff_max-ff_min)
        legend = []
        
        theta = np.deg2rad(np.array(self.data['Theta']))
        phi = np.deg2rad(np.array(self.data['Phi']))
        phi_grid,theta_grid = np.meshgrid(phi, theta)
        
        r = np.reshape(ff_data_renorm,(len(self.data['Theta']),len(self.data['Phi'])))
        
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
        

        plt.show()
        return fig1, ax1
        
    def get_geometry(self,is_antenna_array=True):
        
        time_before = walltime.time()
        print('Exporting Geometry...')
        oEditor = self.aedtapp.odesign.SetActiveEditor("3D Modeler")
        
        #obj is being exported as model units, scaling factor needed for display
        model_units = oEditor.GetModelUnits()
        sf = utils.convert_units(1,'meter',model_units)

        
        bounding_box = oEditor.GetModelBoundingBox()
        self.xmax = float(bounding_box[3])-float(bounding_box[0])
        self.ymax = float(bounding_box[4])-float(bounding_box[1])
        self.zmax = float(bounding_box[5])-float(bounding_box[2])

        
        geo_path = f'{self.aedtapp.results_directory}\\{self.aedtapp.design_name}.results\\geo\\'
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
                for each in self.data['Element_Location']:
                    translated_mesh=cad_mesh.copy()
                    offset_xyz = self.data['Element_Location'][each]*sf
                    if np.abs(2*offset_xyz[0])>self.xmax:#assume array is centere, factor of 2
                        self.xmax=offset_xyz[0]*2
                    if np.abs(2*offset_xyz[1])>self.ymax:#assume array is centere, factor of 2
                        self.ymax=offset_xyz[1]*2
                    translated_mesh.translate(offset_xyz)
                    meshes+=translated_mesh
            else:
                    translated_mesh.translate(self.data['Element_Location'][each]*sf)
                    meshes=translated_mesh
        else:
            print('WARNING: Unable to display CAD Geometry, ' + cad_file + ' is not found')
        
        self.all_max = np.max(np.array([self.xmax,self.ymax,self.zmax]))
        elapsed_time = walltime.time()-time_before
        print(f'Exporting Geometry...Done: {elapsed_time}seconds')
        return meshes

    
    def polar_plot_3d_pyvista(self,ff,
                            qty_str = 'RealizedGain',
                            convert_to_dB=True,
                            cad_mesh=None,
                            position = np.zeros(3),
                            rotation = np.eye(3)):
        
        ff.beamform(phi_scan=0,theta_scan=0)
        ff.get_far_field_mesh(qty_str = qty_str,convert_to_dB=convert_to_dB)
        
            
        engine = MyCustomRoutine2(ff)
        #engine.ff = ff
        
        #plot everything together
        rotation_euler = self.rotationMatrixToEulerAngles(rotation)*180/np.pi


        p = pv.Plotter()

        
        


        p.add_slider_widget(engine.update_phi, rng=[0, 360], value=0, title='Phi', 
                    pointa=(.35, .1), pointb=(.64, .1),style='modern',event_type='always')
        p.add_slider_widget(engine.update_theta, rng=[-180, 180], value=0, title='Theta', 
                    pointa=(.67, .1), pointb=(.98, .1),style='modern',event_type='always')
        
        sargs = dict(height=0.4, vertical=True, position_x=0.05, position_y=0.5)
        #ff_mesh_inst = p.add_mesh(engine.output,smooth_shading=True,cmap="jet",scalar_bar_args=sargs,opacity=0.5)
        #not sure why, but smooth_shading causes this to not update
        ff_mesh_inst = p.add_mesh(engine.output,cmap="jet",scalar_bar_args=sargs)


        if cad_mesh:
            def toggle_vis_ff(flag):
                ff_mesh_inst.SetVisibility(flag)
            def toggle_vis_cad(flag):
                cad.SetVisibility(flag)
            def scale(value=1):
                ff_mesh_inst.SetScale(value,value,value)
                ff_mesh_inst.SetPosition(position)
                ff_mesh_inst.SetOrientation(rotation_euler)
                #p.add_mesh(ff_mesh, smooth_shading=True,cmap="jet")
                return


            ff_toggle = p.add_checkbox_button_widget(toggle_vis_ff, value=True)
            p.add_text('Show Far Fields', position=(70,25), color='black', font_size=12)
            slider_max= int(np.ceil(self.all_max))*2
            scale_slider = p.add_slider_widget(scale, [0, slider_max], title='Scale Plot',value=int(slider_max/2))

            if 'MaterialIds' in cad_mesh.array_names:
                color_display_type = cad_mesh['MaterialIds']
            else:
                color_display_type=None
            cad = p.add_mesh(cad_mesh,scalars=color_display_type,show_scalar_bar=False,opacity=0.5)
    
            cad_toggle = p.add_checkbox_button_widget(toggle_vis_cad, value=True,position=(10,70))
            p.add_text('Show Geometry', position=(70,75), color='black', font_size=12)
        p.show()
        
        return p
        
        
 

    def find_nearest(self,array,value):
        idx = np.searchsorted(array, value, side="left")
        if idx > 0 and (idx == len(array) or math.fabs(value - array[idx-1]) < math.fabs(value - array[idx])):
            return int(array[idx-1])
        else:
            return int(array[idx])

    def rotationMatrixToEulerAngles(self,R) :

        sy = math.sqrt(R[0,0] * R[0,0] +  R[1,0] * R[1,0])
        singular = sy < 1e-6
        if  not singular :
            x = math.atan2(R[2,1] , R[2,2])
            y = math.atan2(-R[2,0], sy)
            z = math.atan2(R[1,0], R[0,0])
        else :
            x = math.atan2(-R[1,2], R[1,1])
            y = math.atan2(-R[2,0], sy)
            z = 0
        return np.array([x, y, z])
    
    def get_new_file_name(self):
        increment=1
        #print("Window size ", p.window_size)
        file_name = self.absolute_path + "\\geo_envelope_overlay" + str(increment) + ".png"
        while os.path.exists(file_name):
            increment+=1
            file_name = self.absolute_path + "\\geo_envelope_overlay" + str(increment) + ".png"
        return file_name
    

class MyCustomRoutine2():
    def __init__(self, ff):
        self.output = ff.mesh
        self._phi=0
        self._theta=0
        # default parameters
        self.ff = ff
        self.qty_str='RealizedGain'
        self.convert_to_dB=True

    def _update_both(self):
        self.ff.beamform(phi_scan=self._phi, theta_scan=self._theta)
        self.ff.get_far_field_mesh(self.qty_str,self.convert_to_dB)
        self.output.overwrite(self.ff.mesh)
        return
    
    def update_phi(self, phi):
        self._phi = phi
        self._update_both()
        
    def update_theta(self, theta):
        self._theta = theta
        self._update_both()
        
