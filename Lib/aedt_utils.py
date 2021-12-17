# -*- coding: utf-8 -*-
"""
Created on Mon Nov 22 11:46:09 2021

interaction between main script and aedt is defined here, assumes
that the desktop is open, uses current active project

@author: asligar
"""


import os


class AEDTutils:
    def __init__(self,aedtapp,project_name='project1',design_name='design1'):
        
        self.aedtapp = aedtapp
        self.project_name = project_name
        self.design_name = design_name
        

        projects = self.aedtapp.project_list
        if self.project_name in projects:
            self.aedtapp.odesktop.SetActiveProject(self.project_name)
        
            designs = self.aedtapp.design_list
            if self.design_name not in designs:
                print(f'Design {self.design_name} not found, please open project/design')
            else:
                self.aedtapp.set_active_design(self.design_name)
        else:
            print(f'Project {self.project_name} not found, please open project')
   
        self.oDesign = self.aedtapp.odesign
        oEditor = self.oDesign.SetActiveEditor("3D Modeler")

        self.model_units = oEditor.GetModelUnits()
        #oEditor.SetModelUnits(["NAME:Units Parameter","Units:=", "meter","Rescale:=", False])
    

        oModule = self.oDesign.GetModule("RadField")
        self.ff_setups_names = oModule.GetSetupNames('Infinite Sphere')


    def release_desktop(self):
        self.aedtapp.release_desktop(close_projects=False, close_on_exit=False)


    

    def add_or_edit_variable(self,name,value):
        self.aedtapp[name]=value
    
