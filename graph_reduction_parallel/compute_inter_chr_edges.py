import argparse
import numpy as np
from graph_reduction_parallel import load_graph
from graph_reduction_parallel import filter_edges
from graph_reduction_parallel import load_patient
from graph_reduction_parallel import compute_nodes_to_merge, merge_nodes, merge_edges
import sys
import networkx as nx
import pandas as pd


def get_args():
    parser = argparse.ArgumentParser('python')

    parser.add_argument('-edge_dir',
                        #default='../../processed_main_graph/final_edge.h5',
                        default='../../test_data/test_2.h5',
                        required=False,
                        help='directory of output edge file.')  

    parser.add_argument('-node_dir',
                        #default='../../processed_main_graph/final_node.csv',
                        default='../../test_data/test_2.csv',
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

    parser.add_argument('-reduced_inter_chr_edges_dir',
                        default='../../chr_data_reduced/',
                        required=False,
                        help='csv file of reduced node table.') 
   
    return parser.parse_args()


def remove_intra_chr_edges(edges_array, nodes_array):
    """
    Remove the intra-chr edges and 
    """
    print(nodes_array)
    print(edges_array)
    print('size of edge table (original): ', edges_array.shape)

    edges_array_copy = edges_array
    for chr in range(1,24):
        print('chr:', chr)
        nodes_array_chr = nodes_array[nodes_array[:,1]==chr, :] # nodes of this chromosome
        nodes_set_chr = set(nodes_array_chr[:,0]) # set of nodes in this chromosome

        edges_to_remove = [] # indices of intra-chr edges to remove
        #print(edges_array_copy)
        for i in range(edges_array_copy.shape[0]):
            row = edges_array_copy[i]
            if row[0] in nodes_set_chr and row[1] in nodes_set_chr: # if it is intra-chr
                edges_to_remove.append(i) # remove
        
        edges_to_remove = np.array(edges_to_remove)
        edges_array_copy = np.delete(edges_array_copy, edges_to_remove, 0) # delete the nodes to merge    

        print('number of intra edges:', edges_to_remove.shape[0])
        print('size of remaining edges table:', edges_array_copy.shape[0])

        print('-----------------------------------------------------')
    return edges_array_copy 


if __name__ == "__main__":
    '''options and input arguments'''
    #np.set_printoptions(threshold=sys.maxsize)
    args = get_args()
    edge_dir = args.edge_dir
    node_dir = args.node_dir
    snps_dir = args.snps_dir
    patient_dir = args.patient_dir
    snp_weight_dir = args.snp_weight_dir
    verbose = int(args.verbose)
    reduced_inter_chr_edges_dir = args.reduced_inter_chr_edges_dir

    print('*****************************************************************')
    print('loading the main graph and preprocessing...')
    nodes_array, edges_array, snp_map, node_id_set = load_graph(node_dir, edge_dir, snps_dir)
    
    '''nodes_array: [node-id, chromosome, chunk_start, chunk_end, has_snp]'''
    nodes_array = load_patient(nodes_array, patient_dir, snp_map, node_id_set, snp_weight_dir, snp_weight_th=0.00016)
    print('total number of nodes in main graph: ', nodes_array.shape[0])
    print('total number of SNP nodes in main graph: ', np.count_nonzero(nodes_array[:,4]==1))

    edges_array = filter_edges(edges_array, th=0.05) # filter out unimportant edges
    print('total number of edges in main graph: ', edges_array.shape[0])

    '''compute inter-chromosome edges to merge'''
    inter_chr_edges = remove_intra_chr_edges(edges_array, nodes_array)
    print('total number of inter-chromosome edges in main graph (before reduction): ', inter_chr_edges.shape[0])

    '''get the dictionary that maps old nodes to new nodes'''
    to_merge = compute_nodes_to_merge(nodes_array)
    nodes_array_reduced, old_to_new_dict = merge_nodes(nodes_array, to_merge, 10)
    inter_chr_edges_reduced = merge_edges(inter_chr_edges, old_to_new_dict, 10)
    print('number of inter_chr_edges: ', inter_chr_edges_reduced.shape[0])

    '''save the inter-chr edge table as a file'''
    patient_code = patient_dir.split('/')[-1]
    patient_code = patient_code.split('.')[0]
    reduced_inter_chr_edges_dir = reduced_inter_chr_edges_dir + '-' + patient_code + '.npy'
    print('saving reduced edge table to: ', reduced_inter_chr_edges_dir)
    np.save(reduced_inter_chr_edges_dir, inter_chr_edges_reduced)  