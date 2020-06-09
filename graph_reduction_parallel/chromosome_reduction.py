"""
Each chromosome is treated as a single graph, and they are reduced seperatedly.
"""
import argparse
import numpy as np
from graph_reduction_parallel import load_graph
from graph_reduction_parallel import filter_edges
from graph_reduction_parallel import load_patient
import sys

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

    parser.add_argument('-snp_weight_dir',
                        default='../../snp_weight/PRS_SNPs',
                        required=False,
                        help='csv file containing weights of snps') 

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

    parser.add_argument('-chr',
                        default=1,
                        required=False,
                        help='which chromosome to compute')       

    return parser.parse_args()


def merge_chr_nodes(nodes_array, chr):
    """
    Merge the nodes of a single chromosome. 
    nodes_array: numpy array containing all nodes of the patient.
    chr: chromosome number.
    """
    '''compute nodes to merge for this chromosome'''
    to_merge = []
    merged_node = []
    for row in nodes_array[nodes_array[:,1]==chr, :]: # for each row in this chromosome
        #print(row)
        if row[4]==0:
            merged_node.append(row[0]) # append node id to merge
        elif row[4]==1:
            if len(merged_node) > 1: # counted as merged on when there are at least 2 nodes in the list
                to_merge.append(merged_node)
            merged_node = []
        else:
            raise ValueError('has_snp column should be either integer of 1 or 0')
    if len(merged_node) > 1: # counted as merged on when there are at least 2 nodes in the list
        to_merge.append(merged_node)

    '''compute merged node table for this chromosome'''
    nodes_array_copy = nodes_array[nodes_array[:,1]==chr, :] # nodes of corresponding chromosomes
    idx_to_remove = []
    for old_nodes in to_merge: # difficult to parallelize because processes will be modifying shared array.
        for old_node in old_nodes:
            idx = np.nonzero(nodes_array_copy[:,0]==old_node)
            idx_to_remove.append(idx[0].item())
    idx_to_remove = np.array(idx_to_remove)
    nodes_array_copy = np.delete(nodes_array_copy, idx_to_remove, 0) # delete the nodes to merge
        
    unchanged_nodes = list(nodes_array_copy[:,0]) # keys also contents for unchanged nodes
    old_to_new_dict_total = dict(zip(unchanged_nodes, unchanged_nodes)) # update the dictionary for edge reduction
    new_to_old_dict_total = dict(zip(unchanged_nodes, unchanged_nodes)) # same for unchanged nodes

    if len(to_merge)>0: 
        old_to_new_dict = {}
        new_nodes = []
        for nodes in to_merge:
            nodes_array_sub = [] # get sub-array of nodes to merge
            for node in nodes: # get the corresponding sub-array node table
                nodes_array_sub.append(nodes_array[np.nonzero(nodes_array[:,0]==node)])
            nodes_array_sub = np.concatenate(nodes_array_sub)
            new_node_id = nodes_array_sub[0, 0] # use first old node id as merged node id
            chromosome =  nodes_array_sub[0, 1] # chromosome of new node
            chunk_start = np.amin(nodes_array_sub[:, 2]) # chunk_start of new node
            chunk_end = np.amax(nodes_array_sub[:, 3])# chunk_end of new node
            for node in nodes: # update the old to new node dictionary
                old_to_new_dict[node] = new_node_id
            new_nodes.append(np.array([new_node_id, chromosome, chunk_start, chunk_end, 0]))
            new_to_old_dict_total[new_node_id] = nodes
            #print(new_to_old_dict_total)
            #print(nodes)
            #print(nodes_array_sub)
            #print('new node id:', new_node_id)
            #print('chr:', chromosome)
            #print('chunk_start:', chunk_start)
            #print('chunk_end:', chunk_end)
        #print('old to new dict:', old_to_new_dict)
        #print('new nodes array:', new_nodes)
        new_nodes = np.vstack(new_nodes)
        old_to_new_dict_total.update(old_to_new_dict)
        nodes_array_copy = np.vstack([nodes_array_copy, new_nodes])
        return nodes_array_copy, old_to_new_dict_total, new_to_old_dict_total
    else:
        return np.empty(shape=(0, 0)), {}, {}


def find_node_set(nodes_array, chr):
    """
    Given a node table and the chromosome number, return a set of node numbers.
    """
    nodes_array = nodes_array[nodes_array[:,1]==chr, :]
    nodes_set = set(nodes_array[:,0])
    return nodes_set


def count_num_inter_chr_nodes(new_nodes, edges_array, new_to_old_dict, old_nodes_set):
    """
    Count number of inter-chromosome nodes for each chromosome.
    """
    num_inter_chr_nodes = 0 # number of inter-chromosome nodes
    for new_node in new_nodes: 
        old_nodes = new_to_old_dict[new_node[0]] # map to corresponding old nodes
        try: # if old_nodes is list
            for old_node in old_nodes:
                idx_source = edges_array[:,0]==old_node # find edges that connected to the original node as source
                idx_target = edges_array[:,1]==old_node # find edges that connected to the original node as target
                idx = idx_source | idx_target    
                idx = np.nonzero(idx)[0] # indexes of rows
                old_edges = edges_array[idx] # array of corresponding old nodes 
                old_edges_source = set(old_edges[:,0]) # sources of old edges
                old_edges_target = set(old_edges[:,1]) # targets of old edges
                old_edges_nodes = set.union(old_edges_source, old_edges_target) # set of nodes (old nodes of this chromosome and other nodes connected with them)
                if not old_edges_nodes.issubset(old_nodes_set): # some node is connected to other chromosome
                    num_inter_chr_nodes += 1
                    break
        except: # if not a list 
            idx_source = edges_array[:,0]==old_nodes # find edges that connected to the original node as source
            idx_target = edges_array[:,1]==old_nodes # find edges that connected to the original node as target
            idx = idx_source | idx_target    
            idx = np.nonzero(idx)[0] # indexes of rows
            old_edges = edges_array[idx] # array of corresponding old nodes 
            old_edges_source = set(old_edges[:,0]) # sources of old edges
            old_edges_target = set(old_edges[:,1]) # targets of old edges
            old_edges_nodes = set.union(old_edges_source, old_edges_target) # set of nodes (old nodes of this chromosome and other nodes connected with them)
            if not old_edges_nodes.issubset(old_nodes_set): # some node is connected to other chromosome
                num_inter_chr_nodes += 1
    return num_inter_chr_nodes


if __name__ == "__main__":
    '''options and input arguments'''
    args = get_args()
    edge_dir = args.edge_dir
    node_dir = args.node_dir
    snps_dir = args.snps_dir
    patient_dir = args.patient_dir
    snp_weight_dir = args.snp_weight_dir
    verbose = int(args.verbose)
    reduced_node_dir = args.reduced_node_dir
    reduced_edge_dir = args.reduced_edge_dir
    reduced_gexf_dir = args.reduced_gexf_dir
    reduced_graph_statistics = args.reduced_graph_statistics
    chr = int(args.chr)

    np.set_printoptions(threshold=sys.maxsize)

    print('computing info for chromosome {}'.format(chr))

    nodes_array, edges_array, snp_map, node_id_set = load_graph(node_dir, edge_dir, snps_dir)
    
    nodes_array = load_patient(nodes_array, patient_dir, snp_map, node_id_set, snp_weight_dir, snp_weight_th=0.00016)
    #print(nodes_array)

    edges_array = filter_edges(edges_array, th=0.05) # filter out unimportant edges
    #print(edges_array)

    print('-----------------------------------------------------------------') 
    print('computing merged nodes of chromosome {}'.format(chr))
    print('number of nodes before merging:', nodes_array[nodes_array[:,1]==chr, :].shape[0])
    old_nodes_set = find_node_set(nodes_array, chr) # set of original nodes of this chromosome        
    #print('set of old nodes: ', old_nodes_set)

    new_nodes, old_to_new_dict, new_to_old_dict = merge_chr_nodes(nodes_array, chr)
    print('number of nodes after merging:', new_nodes.shape[0])
    num_merged_nodes = new_nodes.shape[0]
    print('total number of merged nodes of the patient:', num_merged_nodes)
        
    num_inter_chr_nodes = count_num_inter_chr_nodes(new_nodes, edges_array, new_to_old_dict, old_nodes_set)
    print('number of inter chromosome nodes: ', num_inter_chr_nodes)
    print('-----------------------------------------------------------------')         




