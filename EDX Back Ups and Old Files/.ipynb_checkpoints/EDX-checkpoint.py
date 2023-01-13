import h5py
import matplotlib.pyplot as plt
import numpy as np

def get_hdf_scan(file_dir, file_name):
    fname = file_dir + "\\"+file_name
    # read the array of spectra (in order of bin number)
    with h5py.File(fname, 'r') as f:
        data_spectra = f['dxpMercury:mca1.VAL'][...]
    # read the motor position for the scan of the sample
    # across the beam
    with h5py.File(fname, 'r') as f:
        data_position = f['7bmb1:aero:m2.VAL'][...]
        
    # now you can run a loop over the number of elements in data_position
    # to make a spectrum object at each position
    # and then use that and the position to make
    # all the stuff for the scan object
    
    return(data_spectra)   # this just gives the data array, need
                  # to change is to return a scan object


def bin_to_energy(bins):
    return(bins)
    # need to have everything as np arrays
    
def bin_to_CuKa_angle(bins):
    return(bins)
    # need to have everything as np arrays    
    
    
class spectrum:
    def __init__(self, position, counts):
        self.position = position
        self.bins = np.array(range(0, counts.shape[0] - 1))
        self.counts = np.array(counts)
        self.energies = bin_to_energy(bins)  # need to define this function above
        self.angle = bin_to_CuKa_angle(bins) # need to define this function above
    
    def plot(self):
        plt.plot(self.energies, self. counts)
    
    

class scan:
    def __init__(self, x_coord, spectrum):
        self.x_coords = x_coords
        self.spectrum = spectrum
        
    def plot(self):
        return()
        # make plot wih options for low/hi pixels, and for stack seperation and horizontal perspective
        #plt.plot(self.energies, self. counts)        
    