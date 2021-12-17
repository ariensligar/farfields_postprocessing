# -*- coding: utf-8 -*-
"""
Created on Fri Sep 10 16:02:31 2021

@author: asligar

inculdes some simple functions that are used for various tasks

"""

import numpy as np
import os
import datetime
import json
import csv

def diff(li1, li2): 
    """
    used to return difference between two lists
    commonly used for when HFSS doesn't return the name of objects, for example
    when an stl file is imported, this function can be used to compare list
    of objects before and after import to return the list of imported objects

    returns: difference between lists
    """
    li_dif = [i for i in li1 + li2 if i not in li1 or i not in li2] 
    return li_dif 

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

def round_time(dt=None, roundTo=1):
   """Round a datetime object to any time lapse in seconds
   dt : datetime.datetime object, default now.
   roundTo : Closest number of seconds to round to, default 1 minute.
   """
   if dt == None : dt = datetime.datetime.now()
   seconds = (dt.replace(tzinfo=None) - dt.min).seconds
   rounding = (seconds+roundTo/2) // roundTo * roundTo
   return dt + datetime.timedelta(0,rounding-seconds,-dt.microsecond)
             
def write_csv(data,path='./out.csv',sort_by_pd=False):
    #headings = list(data.keys()) 
    #these heading should always exist
    headings = ['BeamId','Module_Name','PD_Type','EvalSurface','Freq','Averaging_Area','PD_Max','RadiatedPower','Renormalized PD']
    beam_ids = list(data['PD_Max'].keys())
    data['BeamId'] = beam_ids
    #headings.append('BeamId')
    #remove some headings we don't want to output in csv file
    #headings = list(filter(('Paths_To_Raw_Data').__ne__, headings))
    #headings = list(filter(('Paths_To_Avg_Data').__ne__, headings))

    all_rows = []
    if sort_by_pd:
        pd_max_as_list = list(data['PD_Max'].values())
        original_index_vals = list(range(len(beam_ids)))
        zipped = zip(pd_max_as_list,data['BeamId'],original_index_vals)
        sort_zip = list(sorted(zipped,reverse=True))
        for s in sort_zip:
            row_data = []
            for head in headings:
                cell_data= data[head]
                #cell data may be dict, list or single value
                if isinstance(cell_data, list):
                    cell_data=str(cell_data[s[2]])
                elif isinstance(cell_data, dict):
                    cell_data = str(cell_data[s[1]])
                else:
                    cell_data= str(cell_data)
                row_data.append(cell_data)
            all_rows.append(row_data)
    else:
        for n, beam in enumerate(beam_ids):
            row_data = []
            for head in headings:
                cell_data = data[head]
                #cell data may be dict, list or single value
                if isinstance(cell_data, list):
                    cell_data=str(cell_data[n])
                elif isinstance(cell_data, dict):
                    cell_data = str(cell_data[beam])
                else:
                    cell_data= str(cell_data)
                row_data.append(cell_data)
            all_rows.append(row_data)


    with open(path, 'w',newline='') as csvfile: 

        csvwriter = csv.writer(csvfile)  
        csvwriter.writerow(headings) 
        csvwriter.writerows(all_rows)

def dict_with_numpy_to_lists(input_dict):
    '''
    convert a dictionary with numpy arrays to list so they can be dumped to 
    json file. This is tempoarary, works with 2 levels of nested dictionaries.
    

    '''
    for level1 in input_dict.keys():
        if isinstance(input_dict[level1], np.ndarray):
            input_dict[level1] = input_dict[level1].tolist()
        elif isinstance(input_dict[level1], dict):
            for level2 in input_dict[level1].keys():
                if isinstance(input_dict[level1][level2], np.ndarray):
                    input_dict[level1][level2] = input_dict[level1][level2].tolist()
                    
    return input_dict
    





def split_num_units(s):
    if '^' in s:
        s =s.replace('^','')
    s = s.replace(" ","").lower()
    value = s.rstrip('abcdefghijklmnopqrstuvwxyz')
    units = s[len(value):]
    value = float(value)
    return value, units



