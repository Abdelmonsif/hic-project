import argparse
import networkx as nx
import h5py
import numpy as np

def get_args():
    parser = argparse.ArgumentParser('python')

    parser.add_argument('-data_dir',
                        default='../../G_snps_23.gexf',
                        required=False,
                        help='directory of gene expression graph')  

    parser.add_argument('-edge_dir',
                        default='../../G_snps_23.h5',
                        required=False,
                        help='directory of gene expression graph')  
                 
    return parser.parse_args()


class HicGraph:
    def __init__(self, data_dir):
        self.hic_graph = nx.read_gexf(data_dir) 

    def report(self):
        """
        Report info about this graph
        """
        print('number of nodes: ', nx.number_of_nodes(self.hic_graph))
        print('number of edges: ', len(list(self.hic_graph.edges)))
        #print(nx.degree(self.hic_graph))
        #print('node list: ', sorted(list(self.hic_graph.nodes)))

        #print('edge list: ', sorted(list(self.hic_graph.edges)))
        #print('attribute names of each node:', self.node_attribute_list)
        #print('attribute names of each edge:', self.edge_attribute_list)
        print(self.hic_graph.edges['49', '60'])

    def export_edge_list(self, edge_dir):
        """
        Export the edge list along with the attributes of the edges. Each edge 
        has following attributes:
            1. contactCount.
            2. p-value.
            3. q-value.
            4. id.
        """
        

        edge_list = list(self.hic_graph.edges)
        #edge_list = [tuple(map(int, x)) for x in edge_list]
        edge_dict = {} # dictionary that stores edge attributes
        for edge in edge_list:  
            edge_attr = self.hic_graph.edges[edge[0], edge[1]] # retrieve edge attribute   
            edge = tuple(map(int, edge)) # convert edge from string to integer for storage efficiency
            
 

if __name__ == "__main__":
    args = get_args()
    data_dir = args.data_dir
    edge_dir = args.edge_dir

    hic_graph = HicGraph(data_dir)
    hic_graph.report()
    hic_graph.export_edge_list(edge_dir)