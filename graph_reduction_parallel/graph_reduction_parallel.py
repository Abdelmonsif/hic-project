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
    edge_array = edges_df.to_numpy(dtype=float) # convert edge table to numpy array
    return nodes_array, edge_array, snp_map, node_id_set
    

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


def copyto_shared_mem(array, data_type):
    """
    Copy the input numpy array into a shared memory.
    """
    array_shape = array.shape
    print(array_shape)
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

    np.set_printoptions(threshold=sys.maxsize)
    
    '''
    Load the graph as numpy array.
    Node array columns: node_id, chr, chunk_start, chunk_end.
    Edge array columns: source, target, contactCount, p-value, q-value.
    '''
    nodes_array, edges_array, snp_map, node_id_set = load_graph(node_dir, edge_dir, snps_dir)
    #print(nodes_array)
    #print(edges_array)
    
    '''
    Load patient's snp info into nodes_array.
    Node array columns become: node_id, chr, chunk_start, chunk_end, has_snp.
    '''
    nodes_array = load_patient(nodes_array, patient_dir, snp_map, node_id_set)
    #print(nodes_array)

    '''
    Convert original numpy array to shared memory. In the mean time, the original memory 
    is released by re-assignment.
    '''
    nodes_array = copyto_shared_mem(array=nodes_array, data_type='l') # signed long (32-bit)

    '''Compute the nodes to merge'''
    print('computing nodes to merge...')
    to_merge = compute_nodes_to_merge(nodes_array)
    print(to_merge)



