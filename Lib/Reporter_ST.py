# -*- coding: utf-8 -*-
"""
Created on Tue Sep  7 15:50:47 2021

@author: asligar
"""
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from datetime import date
from matplotlib import cm
import pyvista as pv
import math
import copy
import Lib.Utillities as utils

class Report_Module():
    def __init__(self,data):

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
        


            
    def polar_plot_3d_pyvista(self,
                            all_max=1000,
                            qty_str = 'RealizedGain',
                            convert_to_dB=True,
                            cad_mesh=None,
                            position = np.zeros(3),
                            rotation = np.eye(3)):
        if convert_to_dB:
            ff_data = 10*np.log10(self.data[qty_str])
            #renormalize to 0 and 1 
            ff_max_dB = np.max(ff_data)
            ff_min_dB = np.min(ff_data)
            ff_data_renorm = (ff_data-ff_min_dB)/(ff_max_dB-ff_min_dB)
            display_name = f'{qty_str} (dB)'
        else:
            ff_data = self.data[qty_str]
            #renormalize to 0 and 1 
            ff_max = np.max(ff_data)
            ff_min = np.min(ff_data)
            ff_data_renorm = (ff_data-ff_max)/(ff_max-ff_min)
            display_name = f'{qty_str}'
            
        theta = np.deg2rad(np.array(self.data['Theta']))
        phi = np.deg2rad(np.array(self.data['Phi']))
        phi_grid,theta_grid = np.meshgrid(phi, theta)
        
        r_no_renorm = np.reshape(ff_data,(len(self.data['Theta']),len(self.data['Phi'])))
        r = np.reshape(ff_data_renorm,(len(self.data['Theta']),len(self.data['Phi'])))
        
        x = r * np.sin(theta_grid) * np.cos(phi_grid)
        y = r * np.sin(theta_grid) * np.sin(phi_grid)
        z = r * np.cos(theta_grid)
        
        #for color display
        mag = np.ndarray.flatten(r_no_renorm,order='F')
        
        # create a mesh that can be displayed
        ff_mesh = pv.StructuredGrid(x,y,z)
        #ff_mesh.scale(ff_scale)
        #ff_mesh.translate([float(position[0]),float(position[1]),float(position[2])])
        ff_mesh[display_name] = mag
        

        #plot everything together
        rotation_euler = self.rotationMatrixToEulerAngles(rotation)*180/np.pi


        p = pv.Plotter()

        if cad_mesh:
            def toggle_vis_ff(flag):
                ff.SetVisibility(flag)
            def toggle_vis_cad(flag):
                cad.SetVisibility(flag)
            def scale(value=1):
                ff.SetScale(value,value,value)
                ff.SetPosition(position)
                ff.SetOrientation(rotation_euler)
                #p.add_mesh(ff_mesh, smooth_shading=True,cmap="jet")
                return
            
            ff_toggle = p.add_checkbox_button_widget(toggle_vis_ff, value=True)
            ff = p.add_mesh(ff_mesh,smooth_shading=True,cmap="jet")
            slider_max= int(np.ceil(self.all_max))*2
            scale_slider = p.add_slider_widget(scale, [0, slider_max], title='Scale Plot',value=int(slider_max/2))
            
            if 'MaterialIds' in cad_mesh.array_names:
                color_display_type = cad_mesh['MaterialIds']
            else:
                color_display_type=None
            cad = p.add_mesh(cad_mesh,scalars=color_display_type,show_scalar_bar=False,opacity=0.5)
    
            cad_toggle = p.add_checkbox_button_widget(toggle_vis_cad, value=True,position=(10,70))
            

        return p, cad_mesh
        
        
 

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