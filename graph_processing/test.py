import h5py

data_dir = './test.hdf5'
f = h5py.File(data_dir, "w") # create
f = h5py.File(data_dir, "a") # append
edge_list = f.create_group('edge/edgelist')
contactCount = f.create_group('edge/contactCount')
p_values = f.create_group('edge/p_values')
q_values = f.create_group('edge/q_values')
edge_ids = f.create_group('edge/edge_ids')



