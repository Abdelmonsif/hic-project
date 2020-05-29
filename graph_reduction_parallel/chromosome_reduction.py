"""
Each chromosome is treated as a single graph, and they are reduced seperatedly.
"""
import numpy as np
from graph_reduction_parallel import get_args
from graph_reduction_parallel import load_graph
from graph_reduction_parallel import filter_edges
from graph_reduction_parallel import load_patient


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
        return nodes_array_copy
    else:
        return np.empty(shape=(0, 0))

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

    nodes_array, edges_array, snp_map, node_id_set = load_graph(node_dir, edge_dir, snps_dir)
    
    edges_array = filter_edges(edges_array, th=0.05)

    nodes_array = load_patient(nodes_array, patient_dir, snp_map, node_id_set, snp_weight_dir, snp_weight_th=0.00016)
    
    total_merged_nodes = 0
    for chr in range(1,24):
        print('-----------------------------------------------------------------') 
        print('computing merged nodes of chromosome {}'.format(chr))
        print('number of nodes before merging:', nodes_array[nodes_array[:,1]==chr, :].shape[0])
        new_nodes = merge_chr_nodes(nodes_array, chr)
        print('number of nodes after merging:', new_nodes.shape[0])
        total_merged_nodes += new_nodes.shape[0]
    print('-----------------------------------------------------------------')         
    print('total number of merged nodes of the patient:', total_merged_nodes)



