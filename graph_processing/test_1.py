import h5py
import numpy as np

data_dir = './test.hdf5'
f = h5py.File(data_dir, "r") # 
print(f.keys())
edge_grp=f['edge']
test = edge_grp['test']
array_shape = test.shape
array = np.empty([test.shape[0], test.shape[1]])
test.read_direct(array)
print(array)