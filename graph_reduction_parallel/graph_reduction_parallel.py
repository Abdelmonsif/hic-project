"""
Parallel implementation of graph_reduction_np.py
"""
import pandas as pd
import argparse
import json
import h5py
import numpy as np
import sys
from multiprocessing import Pool, RawArray
import time
import math


def get_args():
    parser = argparse.ArgumentParser('python')

    parser.add_argument('-edge_dir',
                        #default='../../processed_main_graph/final_edge.h5',
                        default='../../test_data/test_1.h5',
                        required=False,
                        help='directory of output edge file.')  

    parser.add_argument('-node_dir',
                        #default='../../processed_main_graph/final_node.csv',
                        default='../../test_data/test_1.csv',
                        required=False,
                        help='directory of output edge file.')  

    parser.add_argument('-snps_dir',
                        default='../../snp_map/snp_map.json',
                        required=False,
                        help='location of the snp mapping file.') 

    parser.add_argument('-patient_dir',
                        default='../../patients/BCAC-97446542.csv',
                        required=False,
                        help='location of the patient file.') 

    parser.add_argument('-verbose',
                        default=0,
                        required=False,
                        help='set to 1 for debugging.')

    parser.add_argument('-reduced_node_dir',
                        #default='../../patients_reduced/BCAC-97446542-node.csv',
                        default='../../test_data_reduced/reduced_node_test_1_4.csv',
                        required=False,
                        help='csv file of reduced node table.') 

    parser.add_argument('-reduced_edge_dir',
                        #default='../../patients_reduced/BCAC-97446542-edge.csv',
                        default='../../test_data_reduced/reduced_edge_test_1_4.csv',
                        required=False,
                        help='csv file of reduced node table.') 

    parser.add_argument('-reduced_gexf_dir',
                        #default='../../patients_reduced/BCAC-97446542.gexf',
                        default='../../test_data_reduced/reduced_gexf_test_1_4.csv',
                        required=False,
                        help='csv file of reduced node table.') 

    parser.add_argument('-reduced_graph_statistics',
                        #default='../../patients_reduced/BCAC-97446542-statistics.json',
                        default='../../test_data_reduced/reduced_gexf_statistics_test_1_4.csv',
                        required=False,
                        help='csv file of reduced node table.') 
                 
    return parser.parse_args()


def load_graph(node_dir, edge_dir, snps_dir):
    """
    Load the original main graph. All edge attributes and each node's chromosome, chunk_start and chun_end are loaded.
    The graph reduction algorithm is going to work on the 2 loaded data structures: 'nodes' and 'edge_list'
    """
    '''load the files'''
    with open(snps_dir) as f:
        snp_map = json.load(f)
    edge_list, contactCount, p_values, q_values, edge_ids = load_edge(edge_dir)
    nodes_df = load_node(node_dir) # nodes are pre-sorted according to chr and chunk_start
    
    '''Generate a set of node ids. It is used to intersect with the set of node ids with SNPs.
    Might be deprecated in the future when using the whole main graph.'''
    node_list = nodes_df['node_id'].tolist() # get node list
    node_id_set = set(node_list)

    edge_list = list(map(tuple, edge_list)) # edge list 
    edges_df = pd.DataFrame() # edge table       
    edges_df['source'] = [x[0] for x in edge_list]
    edges_df['target'] = [x[1] for x in edge_list]
    edges_df['contactCount'] = contactCount
    edges_df['p-value'] = p_values
    edges_df['q-value'] = q_values

    nodes_array = nodes_df.to_numpy(dtype=int) # convert node table to numpy array
    edges_array = edges_df.to_numpy(dtype=float) # convert edge table to numpy array
    return nodes_array, edges_array, snp_map, node_id_set
    

def load_edge(edge_dir):
    """
    Load edges from the h5 file, then returns the edges list and edge attributes
    as numpy arrays.
    """
    f = h5py.File(edge_dir, "r")
    edge_list = f['edge/edge_list']
    contactCount = f['edge/contactCount']
    p_values = f['edge/p_values']
    q_values = f['edge/q_values']
    edge_ids = f['edge/edge_ids']
    return edge_list[()], contactCount[()], p_values[()], q_values[()], edge_ids[()] # convert to numpy array


def load_node(node_dir):
    """
    Load nodes from the csv file produced by graph_preprocess.py
    """
    nodes = pd.read_csv(node_dir)
    return nodes


def load_patient(nodes_array, patient_dir, snp_map, node_id_set):
    """
    Load the csv file containing SNPs of a patient, then add the locations of 
    SNPs to the nodes dataframe.
    """
    if nodes_array.shape[1] == 5:
        nodes_array = np.delete(nodes_array, 4, 1) # remove last patient's data if not first patients
    nodes_array = np.insert(nodes_array, 4, 0, axis=1) # add a column to indicate presence of SNPs
    patient_snp = pd.read_csv(patient_dir, sep='	') # load patient SNPs as a dataframe
    snp_cols = [] # list containing all the SNPs of the patient
    snp_cols_1 = patient_snp.columns[(patient_snp == 1).iloc[0]].tolist()
    snp_cols_2 = patient_snp.columns[(patient_snp == 2).iloc[0]].tolist()
    snp_cols.extend(snp_cols_1)
    snp_cols.extend(snp_cols_2)
    snp_locations = [] # find the locations (node ids) of the snps
    num_missing_snp = 0 # number of missing snp for this patient, ignore them
    for snp in snp_cols:
        try:
            snp_locations.append(snp_map[snp][-1]) # last element is node id
        except:
            num_missing_snp += 1
    snp_locations = set(snp_locations) # there are multiple SNPs on a single node, so take the set of this list to remove duplicates
    snp_locations = list(snp_locations.intersection(node_id_set)) # this step may be redundant when using main graph as input
    node_ids = nodes_array[:,0] # first column is node ids
    for snp_location in snp_locations:
        nodes_array[node_ids==snp_location, 4] = 1 # set the rows on 5th 1 if it has snp 
    return nodes_array


def copy_nodes_to_shared_mem(array, data_type):
    """
    Copy the input numpy array into a shared memory.
    """
    array_shape = array.shape
    shared_mem = RawArray(data_type, array_shape[0] * array_shape[1]) # the Raw Array object is 1d
    shared_mem_np = np.frombuffer(shared_mem, dtype=int).reshape(array_shape) # wrap with numpy interface
    np.copyto(shared_mem_np, array)
    return shared_mem_np


def compute_nodes_to_merge(nodes_array):
    """
    Given the node table as a numpy array, compute a list of lists of nodes to merge.
    To do this, we merge the nodes without SNP to nodes with SNP. The result graph will be a new NetworkX 
    graph object. All non-SNP nodes are merged together.
    example output: [[1,2,3],[4,5,6]]
    """
    to_merge = [] # list containing lists of nodes to merge
    for chr in range(1,24): # for each chromosome, can be parallelized
        merged_nodes = []
        merged_node = []
        for row in nodes_array[nodes_array[:,1]==chr, :]: # for each row in this chromosome
            #print(row)
            if row[4]==0:
                merged_node.append(row[0]) # append node id to merge
            elif row[4]==1:
                if len(merged_node) > 1: # counted as merged on when there are at least 2 nodes in the list
                    merged_nodes.append(merged_node)
                merged_node = []
            else:
                raise ValueError('has_snp column should be either integer of 1 or 0')
        if len(merged_node) > 1: # counted as merged on when there are at least 2 nodes in the list
            merged_nodes.append(merged_node)
        to_merge.extend(merged_nodes)
    return to_merge


def compute_nodes_to_merge_parallel(nodes_array, num_processes=4):
    """
    Given the node table as a numpy array, compute a list of lists of nodes to merge.
    To do this, we merge the nodes without SNP to nodes with SNP. The result graph will be a new NetworkX 
    graph object. All non-SNP nodes are merged together.
    example output: [[1,2,3],[4,5,6]]
    """
    chr_list = list(range(1,24))
    with Pool(processes=num_processes, initializer=init_worker, initargs=(nodes_array, nodes_array.shape)) as pool:
        results = pool.map(compute_nodes_to_merge_parallel_worker_func, chr_list)
    to_merge = []
    for l in results:
        to_merge.extend(l)
    return to_merge


def init_worker(X, X_shape):
    """
    Use a global dictionary to create reference to shared dictionary.
    """
    shared_mem_dict['X'] = X
    shared_mem_dict['X_shape'] = X_shape


def compute_nodes_to_merge_parallel_worker_func(chr):
    """
    Each process computes the nodes to merge for a single chromosome.
    """
    nodes_array = np.frombuffer(shared_mem_dict['X'], dtype=int).reshape(shared_mem_dict['X_shape'])
    merged_nodes = []
    merged_node = []
    for row in nodes_array[nodes_array[:,1]==chr, :]: # for each row in this chromosome
        if row[4]==0:
            merged_node.append(row[0]) # append node id to merge
        elif row[4]==1:
            if len(merged_node) > 1: # counted as merged on when there are at least 2 nodes in the list
                merged_nodes.append(merged_node)
            merged_node = []
        else:
            raise ValueError('has_snp column should be either integer of 1 or 0')
    if len(merged_node) > 1: # counted as merged on when there are at least 2 nodes in the list
        merged_nodes.append(merged_node)
    return merged_nodes


def merge_nodes(nodes_array, to_merge, num_processes):
    """
    Single-process implementation to merge nodes. The 'num_process' argument is used to simulate multi-process scenario.
    Return a new node table containing merged nodes.
    nodes_array: input node table.
    to_merge: list of lists of nodes to merge.
    """
    nodes_array_copy = nodes_array # a copy of original node table, we don't want to modify original one
    for old_nodes in to_merge: # difficult to parallelize because processes will be modifying shared array.
        for old_node in old_nodes:
            idx = np.nonzero(nodes_array_copy[:,0]==old_node) # row index
            nodes_array_copy = np.delete(nodes_array_copy, idx, 0) # delete corresponding row

    unchanged_nodes = list(nodes_array_copy[:,0]) # keys also contents for unchanged nodes
    old_to_new_dict = dict(zip(unchanged_nodes, unchanged_nodes)) # update the dictionary for edge reduction
    nodes_per_process = math.ceil(len(to_merge)/num_processes) # size of chunk assigned to each process
    to_merge_chunks = [to_merge[x:x+nodes_per_process] for x in range(0, len(to_merge), nodes_per_process)] # assign each chunk to a process
    for to_merge_chunk in to_merge_chunks: # will be parallelize at this loop level
        #print(to_merge_chunk)
        old_to_new_dict_chunk = {}
        new_nodes_chunk = []
        for nodes in to_merge_chunk:
            nodes_array_sub = [] # get sub-array of nodes to merge from shared memory
            for node in nodes: # get the corresponding sub-array node table
                nodes_array_sub.append(nodes_array[np.nonzero(nodes_array[:,0]==node)])
            nodes_array_sub = np.concatenate(nodes_array_sub)
            new_node_id = nodes_array_sub[0, 0] # use first old node id as merged node id
            chromosome =  nodes_array_sub[0, 1] # chromosome of new node
            chunk_start = np.amin(nodes_array_sub[:, 2]) # chunk_start of new node
            chunk_end = np.amax(nodes_array_sub[:, 3])# chunk_end of new node
            for node in nodes: # update the old to new node dictionary
                old_to_new_dict_chunk[node] = new_node_id
            new_nodes_chunk.append(np.array([new_node_id, chromosome, chunk_start, chunk_end, 0]))
            #print(nodes)
            #print(nodes_array_sub)
            #print('new node id:', new_node_id)
            #print('chr:', chromosome)
            #print('chunk_start:', chunk_start)
            #print('chunk_end:', chunk_end)
        #print('old to new dict:', old_to_new_dict_chunk)
        #print('new nodes array:', new_nodes_chunk)
        new_nodes_chunk = np.vstack(new_nodes_chunk)
        old_to_new_dict.update(old_to_new_dict_chunk)
        nodes_array_copy = np.vstack([nodes_array_copy, new_nodes_chunk])
    return nodes_array_copy, old_to_new_dict


def merge_nodes_parallel(nodes_array, to_merge, num_processes=1):
    """
    Multi-process implementation to merge nodes. The 'num_process' argument is used to simulate multi-process scenario.
    Return a new node table containing merged nodes.
    nodes_array: input node table.
    to_merge: list of lists of nodes to merge.
    """
    '''main process'''
    start_time = time.time()    
    nodes_array_copy = nodes_array # a copy of original node table, we don't want to modify original one
    idx_to_remove = []
    for old_nodes in to_merge: # difficult to parallelize because processes will be modifying shared array.
        for old_node in old_nodes:
            idx = np.nonzero(nodes_array_copy[:,0]==old_node)
            idx_to_remove.append(idx[0].item())
    idx_to_remove = np.array(idx_to_remove)
    nodes_array_copy = np.delete(nodes_array_copy, idx_to_remove, 0)
    unchanged_nodes = list(nodes_array_copy[:,0]) # keys also contents for unchanged nodes
    old_to_new_dict = dict(zip(unchanged_nodes, unchanged_nodes)) # update the dictionary for edge reduction
    nodes_per_process = math.ceil(len(to_merge)/num_processes) # size of chunk assigned to each process
    to_merge_chunks = [to_merge[x:x+nodes_per_process] for x in range(0, len(to_merge), nodes_per_process)] # assign each chunk to a process
    print('time of main process:')
    report_elapsed_time(start_time)  
    
    '''child processes'''
    with Pool(processes=num_processes, initializer=init_worker, initargs=(nodes_array, nodes_array.shape)) as pool:
        results = pool.map(merge_nodes_parallel_worker_func, to_merge_chunks)

    '''join to main process'''
    for result in results:
        old_to_new_dict.update(result[0])
        nodes_array_copy = np.vstack([nodes_array_copy, result[1]])

    return nodes_array_copy, old_to_new_dict   


def merge_nodes_parallel_worker_func(to_merge_chunk):
    """
    Merge the nodes in given chunk and return a sub-node table for reduced nodes.
    """
    nodes_array = np.frombuffer(shared_mem_dict['X'], dtype=int).reshape(shared_mem_dict['X_shape']) # declare shared memory
    old_to_new_dict_chunk = {}
    new_nodes_chunk = []
    for nodes in to_merge_chunk:
        nodes_array_sub = [] # get sub-array of nodes to merge from shared memory
        for node in nodes: # get the corresponding sub-array node table
            nodes_array_sub.append(nodes_array[np.nonzero(nodes_array[:,0]==node)])
        nodes_array_sub = np.concatenate(nodes_array_sub)
        new_node_id = nodes_array_sub[0, 0] # use first old node id as merged node id
        chromosome =  nodes_array_sub[0, 1] # chromosome of new node
        chunk_start = np.amin(nodes_array_sub[:, 2]) # chunk_start of new node
        chunk_end = np.amax(nodes_array_sub[:, 3])# chunk_end of new node
        for node in nodes: # update the old to new node dictionary
            old_to_new_dict_chunk[node] = new_node_id
        new_nodes_chunk.append(np.array([new_node_id, chromosome, chunk_start, chunk_end, 0]))
        #print(nodes)
        #print(nodes_array_sub)
        #print('new node id:', new_node_id)
        #print('chr:', chromosome)
        #print('chunk_start:', chunk_start)
        #print('chunk_end:', chunk_end)
    #print('old to new dict:', old_to_new_dict_chunk)
    #print('new nodes array:', new_nodes_chunk)
    new_nodes_chunk = np.vstack(new_nodes_chunk)
    return (old_to_new_dict_chunk, new_nodes_chunk) # can only return one object


def copy_edges_to_shared_mem(array, data_type):
    """
    Copy the input numpy array into a shared memory.
    """
    array_shape = array.shape
    shared_mem = RawArray(data_type, array_shape[0] * array_shape[1]) # the Raw Array object is 1d
    shared_mem_np = np.frombuffer(shared_mem, dtype=np.single).reshape(array_shape) # wrap with numpy interface
    np.copyto(shared_mem_np, array)
    return shared_mem_np


def merge_edges(edges_array, old_to_new_dict, num_processes):
    """
    Merge the edge table according to merged nodes.
    Columns of edge table: source, target, contactCount, p-value, q-value.
    """
    edges_array_reduced = edges_array # make a copy
    node_map = np.vectorize(old_to_new_dict.get) # mapping function to map from old node ids to new node ids
    new_source = node_map(edges_array_reduced[:,0]) # convert source from old to new
    new_target = node_map(edges_array_reduced[:,1]) # convert target from old to new
    
    source_target_mat = np.hstack((new_source.reshape(len(new_source), 1), new_target.reshape(len(new_target), 1))) # two columns of source-target pairs
    source_target_mat.sort(axis=1)# make source always <= target

    edges_array_reduced[:,0] = source_target_mat[:, 0] # put new sources in edge table
    edges_array_reduced[:,1] = source_target_mat[:, 1] # put new targets in edge table
    
    # remove self-loops
    idx_to_remove =np.nonzero(edges_array_reduced[:,0]==edges_array_reduced[:,1])[0]
    edges_array_reduced = np.delete(edges_array_reduced, idx_to_remove, 0)

    source_target_pairs = list(set([tuple(x) for x in edges_array_reduced[:,0:2]])) # set of source-target pairs (source always <= target)
    new_edges = [] # list of new edges
    for source_target_pair in source_target_pairs: # parallelize this loop
        idx_source = edges_array_reduced[:,0]==source_target_pair[0]
        idx_target = edges_array_reduced[:,1]==source_target_pair[1]
        pair_idx = idx_source & idx_target
        pair_idx = np.nonzero(pair_idx)[0] # indexes of rows that has source and target in this pair
        old_edges = edges_array_reduced[pair_idx] # array of corresponding old nodes 
        new_edge = np.median(old_edges, axis=0)
        new_edges.append(new_edge)
    edges_array_reduced = np.vstack(new_edges) # put the new edges in one numpy array
    print(edges_array_reduced.shape)
    return edges_array_reduced


def report_elapsed_time(start):
    end = time.time()   
    time_elapsed = end - start
    hour = time_elapsed//3600
    time_elapsed = time_elapsed - hour * 3600
    minute = time_elapsed//60
    time_elapsed = time_elapsed - minute * 60
    print('{}:{}:{}'.format(int(hour), int(minute), round(time_elapsed)))


shared_mem_dict={} # dictionary pointing to shared memory
if __name__ == "__main__":
    '''options and input arguments'''
    args = get_args()
    edge_dir = args.edge_dir
    node_dir = args.node_dir
    snps_dir = args.snps_dir
    patient_dir = args.patient_dir
    verbose = int(args.verbose)
    reduced_node_dir = args.reduced_node_dir
    reduced_edge_dir = args.reduced_edge_dir
    reduced_gexf_dir = args.reduced_gexf_dir
    reduced_graph_statistics = args.reduced_graph_statistics

    #np.set_printoptions(threshold=sys.maxsize)
    
    '''
    Load the graph as numpy array.
    Node array columns: node_id, chr, chunk_start, chunk_end.
    Edge array columns: source, target, contactCount, p-value, q-value.
    '''
    print('/**************************************************************/')
    print('Loading the graph......')
    start_time = time.time()
    nodes_array, edges_array, snp_map, node_id_set = load_graph(node_dir, edge_dir, snps_dir)
    report_elapsed_time(start_time)   
    #print(nodes_array)
    #print(edges_array)
    
    '''
    Load patient's snp info into nodes_array.
    Node array columns become: node_id, chr, chunk_start, chunk_end, has_snp.
    '''
    print('/**************************************************************/')
    print('Loading patient info and add to graph......')
    start_time = time.time()
    nodes_array = load_patient(nodes_array, patient_dir, snp_map, node_id_set)
    report_elapsed_time(start_time)    

    '''
    Convert original numpy array to shared memory. In the mean time, the original memory 
    is released by re-assignment.
    '''    
    print('/**************************************************************/')
    print('Copying nodes_array to shared memory......')
    start_time = time.time()    
    nodes_array = copy_nodes_to_shared_mem(array=nodes_array, data_type='l') # signed long (32-bit)
    report_elapsed_time(start_time)   


    '''Compute the nodes to merge'''
    print('/**************************************************************/')
    print('Computing nodes to merge with single process......')
    start_time = time.time()    
    to_merge = compute_nodes_to_merge(nodes_array)
    report_elapsed_time(start_time)   
    
    print('/**************************************************************/')
    print('Computing nodes to merge with multiple process......')
    start_time = time.time()    
    to_merge_parallel = compute_nodes_to_merge_parallel(nodes_array=nodes_array, num_processes=20)
    report_elapsed_time(start_time)   
    
    print('testing if single process version matches multi process version:')
    print(to_merge.sort()==to_merge_parallel.sort())

    print('/**************************************************************/')
    print('Generating merged node table (single process)......')
    start_time = time.time() 
    nodes_array_reduced, old_to_new_dict = merge_nodes(nodes_array, to_merge, 10)
    report_elapsed_time(start_time)  

    print('/**************************************************************/')
    print('Generating merged node table (multiple process)......')
    start_time = time.time() 
    nodes_array_reduced_parallel, old_to_new_dict_parallel = merge_nodes_parallel(nodes_array, to_merge, 10)
    report_elapsed_time(start_time)  

    print('testing if single process version matches multi process version (reduced node table)')
    nodes_array_reduced = nodes_array_reduced[nodes_array_reduced[:,0].argsort()]
    nodes_array_reduced_parallel = nodes_array_reduced_parallel[nodes_array_reduced_parallel[:,0].argsort()]
    print(np.array_equal(nodes_array_reduced, nodes_array_reduced_parallel))

    print('testing if single process version matches multi process version (old to new node dictionary)')
    print(old_to_new_dict == old_to_new_dict_parallel)

    '''
    Convert original numpy array to shared memory. In the mean time, the original memory 
    is released by re-assignment.
    '''    
    print('/**************************************************************/')
    print('Copying edges_array to shared memory......')
    start_time = time.time()    
    edges_array = copy_edges_to_shared_mem(array=edges_array, data_type='f') # 32-bit float
    report_elapsed_time(start_time)

    print('/**************************************************************/')
    print('Generating merged edge table (single process)......')
    start_time = time.time() 
    edges_array_reduced = merge_edges(edges_array, old_to_new_dict, 10)
    report_elapsed_time(start_time)  
