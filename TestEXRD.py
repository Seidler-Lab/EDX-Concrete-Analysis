import h5py
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import seaborn as sns
from scipy.misc import derivative
from pathlib import Path
from matplotlib.offsetbox import AnchoredText
from mpl_toolkits import mplot3d

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
        hdf_file = h5py.File(directory + files[y], 'r')
        for z in range(len(raw_df.index)):
            if raw_df['ScanNum'].values[z] == hdf_file.attrs['Scan Number']:
                raw_df['FileName'].values[z] = directory+files[y]
    
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
    
    
directory = '/Users/anthonygironda/Documents/MDA_Scans_Seidler/'
meta = '7bmb1_10.19.22_EDXindex.csv'

sample_df = gen_spec_df(directory, meta)
#sample_df.head()

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
            
scan_n = 1106

s1 = scan_instant_helper(sample_df, scan_n, 13, False)
s2 = scan_instant_helper(sample_df, scan_n, 12.75, False)
s3 = scan_instant_helper(sample_df, scan_n, 12.5, False)
s4 = scan_instant_helper(sample_df, scan_n, 12.25, False)
s5 = scan_instant_helper(sample_df, scan_n, 12, False)
s6 = scan_instant_helper(sample_df, scan_n, 11.75, False)
s7 = scan_instant_helper(sample_df, scan_n, 11.5, False)
s8 = scan_instant_helper(sample_df, scan_n, 11.25, False)
s9 = scan_instant_helper(sample_df, scan_n, 11, False)
s10 = scan_instant_helper(sample_df, scan_n, 10.75, False)
s11 = scan_instant_helper(sample_df, scan_n, 10.5, False)
s12 = scan_instant_helper(sample_df, scan_n, 10.25, False)
s13 = scan_instant_helper(sample_df, scan_n, 10, False)
s14 = scan_instant_helper(sample_df, scan_n, 9.75, False)
s15 = scan_instant_helper(sample_df, scan_n, 9.5, False)
s16 = scan_instant_helper(sample_df, scan_n, 9.25, False)
s17 = scan_instant_helper(sample_df, scan_n, 9, False)
s18 = scan_instant_helper(sample_df, scan_n, 8.75, False)
s19 = scan_instant_helper(sample_df, scan_n, 8.5, False)
s20 = scan_instant_helper(sample_df, scan_n, 8.25, False)
s21 = scan_instant_helper(sample_df, scan_n, 8, False)
s22 = scan_instant_helper(sample_df, scan_n, 7.75, False)
s23 = scan_instant_helper(sample_df, scan_n, 7.5, False)
s24 = scan_instant_helper(sample_df, scan_n, 7.25, False)
s25 = scan_instant_helper(sample_df, scan_n, 7, False)
s26 = scan_instant_helper(sample_df, scan_n, 6.75, False)
s27 = scan_instant_helper(sample_df, scan_n, 6.5, False)
s28 = scan_instant_helper(sample_df, scan_n, 6.25, False)
s29 = scan_instant_helper(sample_df, scan_n, 6, False)

slist = [s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12,s13,s14,s15,s16,s17,s18,s19,s20,s21,s22,s23,s24,s25,s26,s27,s28,s29]

class sample:
    
    recipe_ph = None
    recipe_steel = None
    
    aging_temp = None
    aging_atm = None
    aging_humidity = None
    aging_time = None
    
    def __init__(self, sample_name, scan_set):
        self.name = sample_name
        self.scans = scan_set
        self.scan_len = len(self.scans)
            
    #Set attributes of the class, the instantiation helper function should do most of this
    def set_ph(self, hi_low_string): self.recipe_ph = hi_low_string
    def set_steel_fiber(self, true_false): self.recipe_steel = true_false
    def set_temp(self, temperature): self.aging_temp = temperature
    def set_atm(self, atmosphere): self.aging_atm = atmosphere
    def set_humidity(self, rh): self.aging_humidity = rh
    def set_time(self, days): self.aging_time = days

    #Return values of each attribute
    def get_ph(self): return (self.recipe_ph)
    def get_steel_fiber(self): return (self.recipe_steel)    
    def get_temp(self): return (self.aging_temp)    
    def get_atm(self): return (self.aging_atm)
    def get_humidity(self): return (self.aging_humidity)
    def get_time(self): return (self.aging_time)
    def get_len(self): return self.scan_len 

    def get_positions(self):
        scan_pos = []
        for a in range(self.scan_len):
            if (self.scans[a]).beam_position == None:
                scan_pos.append(0)
            else:
                scan_pos.append((self.scans[a]).beam_position)
        return scan_pos
      
    #Quick plotting functions of the class:
    def quick_plot(self):
        length = self.get_len()
        x_data_bins = list(range(len((self.scans[0]).spectra)))
        fig, ax = plt.subplots()#(figsize=(5, 2.7), layout='constrained')
        for g in range(length):
            ax.plot(x_data_bins, (self.scans[g]).spectra, linewidth=0.5)
        ax.grid(True)
        return (fig, ax)
    
    def quick_wire(self):
        length = self.get_len()
        x_data_bins = list(range(len((self.scans[0]).spectra)))
        pos_array = self.get_positions()
        dataset = []
        for r in range(length):
            tup = self.scans[r].spectra, pos_array[r]
            dataset.append(tup)
        return dataset
    
    def quick_stack(self):
        return
        
    def plot(self):
        return
    
    def stackplot(self):
        return
    
    def plot_d(self):
        return
    
    def stackplot_d(self):
        return
    
hpbl2514 = sample('High pH BL 25C 14d', slist)
hpbl2514.set_ph('High')
hpbl2514.set_steel_fiber(False)
hpbl2514.set_temp(25)
hpbl2514.set_atm('100% CO2')
hpbl2514.set_time(14)

hell = hpbl2514.quick_wire()
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as plt

fig = plt.figure()
ax = fig.add_subplot(projection= '3d')

# def plot2din3d(x,y,z):
#     ax.plot(x, y, zs=z, zdir='z')
#     a2d_col_obj = ax.fill_between(x, 0.5, y, step='pre', alpha=0.1) 
#     ax.add_collection3d(a2d_col_obj, zs = z, zdir = 'z')
    
# Plot a sin curve using the x and y axes.
for z in range(len(hell)):
    x = list(range(len(hell[z][0])))
    y = hell[z][0]
    zed = hell[z][1]
    ax.plot(x, y, zed, linewidth = 0.3, color = 'tab:blue')
    
#ax.legend()
ax.set_xlim(750, 1750)
ax.set_ylim(0, 500)
ax.set_zlim(5, 14)
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')

ax.view_init(150, -90)

plt.show()