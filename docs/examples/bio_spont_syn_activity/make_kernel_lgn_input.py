from os.path import join
import numpy as np
import h5py

datafolder = "spikes_inputs"
input_file_name = "lgn_spikes_kernel.h5"

kernel_spike_time = 100.0

f = h5py.File(join(datafolder, input_file_name), 'a')

f["spikes"]["lgn"]["timestamps"][:] = kernel_spike_time

f.close()