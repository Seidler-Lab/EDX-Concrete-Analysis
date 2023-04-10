
import h5py
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import seaborn as sns
from scipy.misc import derivative
from pathlib import Path


class scan:
    
    normalization = 1
    cal_feat_1 = (275, 2292.1)
    cal_feat_2 = (934, 8397.6)
    units = 'bins'
    
    def __init__(self, spectra, position_from_center):
        self.spectra = spectra
        self.beam_position = position_from_center

    def set_energy_scale(self, tuple1, tuple2):
        self.cal_feat_1 = tuple1
        self.cal_feat_2 = tuple2
        
    def set_units(self, unit):
        if unit == 'energy' or unit == 'bins' or unit == '2 theta':
            self.units = unit
        else: print('Unit must be set to "energy" , "2 theta", or "bins"')

    def get_energy_scale(self): return (self.cal_feat_1, self.cal_feat_2)

    def get_units(self): return (self.units)
         
    def quick_plot(self):
        spec_data = self.spectra
        x_dat = list(range(len(spec_data)))
        fig, ax = plt.subplots()
        ax.plot(x_dat, spec_data, linewidth=2.0)
        ax.grid(True);
        return (fig, ax)
    
    def quick_scatter(self):
        spec_data = self.spectra
        x_dat = list(range(len(spec_data)))
        fig, ax = plt.subplots()
        ax.scatter(x_dat, spec_data)
        ax.grid(True);
        return (fig, ax)
    
    def quick_scatter1(self):
        spec_data = self.spectra
        x_dat = list(range(len(spec_data)))
        x_dat_en = []
        x1 = self.cal_feat_1[0]
        y1 = self.cal_feat_1[1]
        x2 = self.cal_feat_2[0]
        y2 = self.cal_feat_2[1]
        for a in range(len(x_dat)):
            lin_interp = (y2-y1)*(a-x1)/(x2-x1) + y1
            x_dat_en.append(lin_interp)
        fig, ax = plt.subplots(figsize = (10,6))
        ax.plot(x_dat_en, spec_data)
        ax.grid(True);
        return (fig, ax)
    
    def bins_to_energy(self, x_data_array):
        x_data_en = []
        x1 = self.cal_feat_1[0]
        y1 = self.cal_feat_1[1]
        x2 = self.cal_feat_2[0]
        y2 = self.cal_feat_2[1]
        for a in range(len(x_data_array)):
            lin_interp = (y2-y1)*(x_data_array[a]-x1)/(x2-x1) + y1
            x_data_en.append(lin_interp)
        return x_data_en

    def bins_to_energy2(self):
        spec_data = self.spectra
        x_data = list(range(len(spec_data)))
        x_dat_en = []
        x1 = self.cal_feat_1[0]
        y1 = self.cal_feat_1[1]
        x2 = self.cal_feat_2[0]
        y2 = self.cal_feat_2[1]
        for a in range(len(x_data)):
            lin_interp = (y2-y1)*(a-x1)/(x2-x1) + y1
            x_dat_en.append(lin_interp)
        return x_data_en
    
    
#Make a dataframe from the sample directory. Good to visualize data. You should 
#feed this dataframe to another function that will make an array of objects
def gen_spec_df(folder_dir, meta_dic):
    
    #initialize a df of the CSV sample/scan info.
    #This should probably be something I iron out with path (it'd be a '\' on PC) but it's fine for now
    raw_df = pd.read_csv(folder_dir + "/" + meta_dic)
    
    #get the hdf filename for every file in the directory
    files = []
    for x in os.listdir(folder_dir):
        if x.endswith(".hdf5"):
            files.append(x)
    
    #Add a column for filename in meta_df that is the same length of the meta_df
    edxname = [''] * len(raw_df.index)
    raw_df['FileName'] = edxname
    
    #Find each file in the meta_df by comparing scan number. Add file name 
    for y in range(len(files)):
        hdf_file = h5py.File(folder_dir + files[y], 'r')
        for z in range(len(raw_df.index)):
            if raw_df['ScanNum'].values[z] == hdf_file.attrs['Scan Number']:
                raw_df['FileName'].values[z] = folder_dir+files[y]
    
    #Add a column for EDX data arrays
    edxd = [''] * len(raw_df.index)
    raw_df['Spectra'] = edxd
    raw_df['Position'] = edxd
    
    #Scan through filenames and add position and spectra data arrays to df
    #If there's none, set those values equal to some string or warning
    for h in range(len(raw_df.index)):
        datahdf = h5py.File(raw_df['FileName'].values[h], 'r')
        positionflag = 'No Position Data (uh oh)'
        spectraflag = 'No Spectra Data (uh oh)'
        with h5py.File(raw_df['FileName'].values[h], 'r') as f:
            for x in range(len(f.keys())):
                if list(f.keys())[x] == '7bmb1:aero:m2.VAL':
                    raw_df['Position'].values[h] = f['7bmb1:aero:m2.VAL'][...]
                    positionflag = 'yep'
                if list(f.keys())[x] == 'dxpMercury:mca1.VAL':
                    raw_df['Spectra'].values[h] = f['dxpMercury:mca1.VAL'][...]
                    spectraflag = 'yep'
        if positionflag != 'yep':
            raw_df['Position'].values[h] = []
        if spectraflag != 'yep':
            raw_df['Spectra'].values[h] = []

    return raw_df


#When bool_pos_given_as_index is True, the pos instantiation variable is an index. Otherwise, when false,
#the pos instantiation variable is an actual position value, like -13.0 or 6.75. It also catches
#an edge case when the position array is empty, such as the calibration scan 1069. Things may behave funny when
#position arrays end in duplicates of 0s or any duplicates in general

#There's potentially a problem for if the spectra passed to the scan is a list or a numpy array, be aware.
def scan_instant_helper(df, scannum, pos, bool_pos_given_as_index):
    errorflag = True
    pos_i = None
    scan_index = (df[df['ScanNum']==scannum].index.values)[0]
    pos_array = df['Position'][scan_index]
    
    #edge case checker, if the position array of a given scan number is empty
    if ((pos == 0) and len(pos_array)==0):
        spectra_set = df['Spectra'][scan_index]
        spectra = spectra_set
        return scan(spectra, None)

    #If the bool is given as true, treat the pos variable like an index and not a value
    #Check first that the pos as an index is within the range of the scan's position array
    if bool_pos_given_as_index == True:
        if pos > (len(pos_array)-1): 
            return print('ERROR: Index is out of range of position array')
        else:
            pos_i = pos
            pos_value = abs(df['Position'][scan_index][pos_i])
            spectra_set = df['Spectra'][scan_index]
            spectra = spectra_set[pos_i]
            return scan(spectra, pos_value)
    
    #Otherwise, the bool was false and the pos variable was given as an actual position in the array
    #Check first that the pos value is actually in the array, and if so instantiate the scan with appropriate info
    else:
        for a in range(len(pos_array)):
            pos_val = pos_array[a]
            #This old logic in the 'or' statement was to catch the empty edge case, it should be unnecessary now
            if abs(pos_array[a])==abs(pos) or ((pos == 0) and len(pos_array)==0): 
                errorflag = False
                pos_i = a
                
        if errorflag == True: print('ERROR: Position is out of range')
        else:
            spectra_set = df['Spectra'][scan_index]
            spectra = spectra_set[pos_i]
            position_val = abs(df['Position'][scan_index][pos_i])
            return scan(spectra, position_val)