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
                        default='../../G_snps_23.hdf5',
                        required=False,
                        help='directory of gene expression graph')  
                 
    return parser.parse_args()


class HicGraph:
    def __init__(self, data_dir):
        self.hic_graph = nx.read_gexf(data_dir) 
        self.num_nodes = nx.number_of_nodes(self.hic_graph)
        self.num_edges = len(list(self.hic_graph.edges))
        
    def report(self):
        """
        Report info about this graph
        """
        print('number of nodes: ', self.num_nodes)
        print('number of edges: ', self.num_edges)

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
        f = h5py.File(edge_dir, "w") # create
        f = h5py.File(edge_dir, "a") # append
        edge_grp = f.create_group('edge')
        
        edge_list_array = np.empty([self.num_edges, 2]) # allocate memory for edge list (num_edges*2)
        contactCount_array = np.empty(self.num_edges) # allocate memory for contactCount (num_edges)
        p_value_array = np.empty(self.num_edges) # allocate memory for p-values (num_edges)
        q_value_array = np.empty(self.num_edges) # allocate memory for q-values (num_edges)
        edge_id_array = np.empty(self.num_edges) # allocate memory for edge ids (num_edges)

        edge_list = list(self.hic_graph.edges)
        i = 0
        for edge in edge_list: # this loop can be parallelized in the future  
            edge_attr = self.hic_graph.edges[edge[0], edge[1]] # retrieve edge attribute   
            edge = tuple(map(int, edge)) # convert edge from string to integer for storage efficiency
            print(edge)
            print(edge_attr)
            edge_list_array[i] = np.array(edge) # put edge pair in a row of edge list array
            contactCount_array[i] = edge_attr['contactCount'] # put contactCount in ith element of contactCount_array
            p_value_array[i] = edge_attr['p-value'] # put p-value in the ith element of p_value_array
            q_value_array[i] = edge_attr['q-value'] # put q-value in the ith element of q_value_array
            edge_id_array[i] = edge_attr['id']# put edge id in the ith element of edge_id_array
            i += 1
        #-------------------- end of for loop ---------------------
        print(edge_list_array)
        print(contactCount_array)
        print(p_value_array)
        print(q_value_array)
        print(edge_id_array)
        edge_grp.create_dataset('edge_list', data=edge_list_array) # write to hdf5    
        edge_grp.create_dataset('contactCount', data=contactCount_array)
        edge_grp.create_dataset('p_values', data=p_value_array)
        edge_grp.create_dataset('q_values', data=q_value_array)
        edge_grp.create_dataset('edge_ids', data=edge_id_array)

if __name__ == "__main__":
    args = get_args()
    data_dir = args.data_dir
    edge_dir = args.edge_dir

    hic_graph = HicGraph(data_dir)
    hic_graph.report()
    hic_graph.export_edge_list(edge_dir)