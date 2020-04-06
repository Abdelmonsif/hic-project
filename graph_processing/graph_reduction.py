"""
Reduce the size of the graph by contracting edges.
"""
import pandas as pd
import argparse
from hic_graph import HicGraph

def get_args():
    parser = argparse.ArgumentParser('python')

    parser.add_argument('-edge_dir',
                        #default='../../test_data/G_snps_23_edge.h5',
                        default='../../test_data/test_1.h5',
                        required=False,
                        help='directory of output edge file.')  

    parser.add_argument('-node_dir',
                        #default='../../test_data/G_snps_23_node.csv',
                        default='../../test_data/test_1.csv',
                        required=False,
                        help='directory of output edge file.')  

    parser.add_argument('-snps_dir',
                        default='../../snp_map/snp_map.json',
                        required=False,
                        help='location of the snp mapping file.') 

    parser.add_argument('-patient_dir',
                        default='../../patient/BCAC-97446542.csv',
                        required=False,
                        help='location of the patient file.') 

    parser.add_argument('-verbose',
                        default=0,
                        required=False,
                        help='set to 1 for debugging.')

    parser.add_argument('-reduced_node_dir',
                        default='../../test_data_reduced/reduced_node_test_1.csv',
                        required=False,
                        help='csv file of reduced node table.') 

    parser.add_argument('-reduced_edge_dir',
                        default='../../test_data_reduced/reduced_edge_test_1.csv',
                        required=False,
                        help='csv file of reduced node table.') 
                 
    return parser.parse_args()


if __name__ == "__main__":
    pd.set_option("display.max_columns", 8)
    pd.set_option('display.width', 200)
    args = get_args()
    edge_dir = args.edge_dir
    node_dir = args.node_dir
    snps_dir = args.snps_dir
    patient_dir = args.patient_dir
    verbose = int(args.verbose)
    reduced_node_dir = args.reduced_node_dir
    reduced_edge_dir = args.reduced_edge_dir
    hic_graph = HicGraph(edge_dir, node_dir, snps_dir, None, verbose)
    hic_graph.load_graph()
    hic_graph.load_patient(patient_dir)
    hic_graph.graph_reduce_1()
    hic_graph.export_reduced_graph(reduced_node_dir, reduced_edge_dir)




