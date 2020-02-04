"""
Pre-process the gexf file for hi-c graph for storage efficiency. There are 2 output files, 
one for edges and one for the nodes.
"""
import argparse
import networkx as nx
import h5py
import numpy as np
import pandas as pd
import time 

def get_args():
    parser = argparse.ArgumentParser('python')

    parser.add_argument('-data_dir',
                        default='../../test_data/test.gexf',
                        required=False,
                        help='directory of original hi-c graph')  

    parser.add_argument('-edge_dir',
                        default='../../test_data/test.h5',
                        required=False,
                        help='directory of output edge file.')  

    parser.add_argument('-node_dir',
                        default='../../test_data/test.csv',
                        required=False,
                        help='directory of output edge file.')  

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
        
        #print(self.hic_graph.edges['49', '60'])

    def export_graph(self, edge_dir, node_dir):
        """
        Export the information of the graph to an h5 file.
        """

        self.__export_edge_list(edge_dir)
        self.__export_node_list(node_dir)

    def __export_edge_list(self, out_dir):
        """
        Export the edge list along with the attributes of the edges to an h5 file. Each edge 
        has following attributes:
            1. contactCount.
            2. p-value.
            3. q-value.
            4. id.
        """
        start=time.time()
        f = h5py.File(out_dir, "w") # create
        f = h5py.File(out_dir, "a") # append
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
            edge_list_array[i] = np.array(edge) # put edge pair in a row of edge list array
            contactCount_array[i] = edge_attr['contactCount'] # put contactCount in ith element of contactCount_array
            p_value_array[i] = edge_attr['p-value'] # put p-value in the ith element of p_value_array
            q_value_array[i] = edge_attr['q-value'] # put q-value in the ith element of q_value_array
            edge_id_array[i] = edge_attr['id']# put edge id in the ith element of edge_id_array
            i += 1
        #-------------------- end of for loop ---------------------
        edge_grp.create_dataset('edge_list', data=edge_list_array, dtype='i') # write to hdf5    
        edge_grp.create_dataset('contactCount', data=contactCount_array, dtype='i')
        edge_grp.create_dataset('p_values', data=p_value_array, dtype='f')
        edge_grp.create_dataset('q_values', data=q_value_array, dtype='f')
        edge_grp.create_dataset('edge_ids', data=edge_id_array, dtype='i')

        time_elapsed = time.time()-start 
        print('time elapsed for exporting edges: ', time_elapsed)

    def __export_node_list(self, out_dir):
        """
        Export the info of nodes to an h5 file. Each line of the csv file has the following attributes:
            1. node id (integer starts with 1).
            2. chromosome number (integer from 1 to 23).
            3. chunk-start (start position of the chunk represented by this node, integer).
            4. chunk-end (end position of the chunk represented by this node, integer).
        """
        start = time.time()
        index = range(self.num_nodes)
        columns = ['node_id', 'chr', 'chunk_start', 'chunk_end']
        df = pd.DataFrame(index=index, columns=columns)
        node_list = list(self.hic_graph.nodes)
        i = 0
        for node in node_list:
            node_attr = self.hic_graph.nodes[node]
            assert(node_attr['label']==node) # node id and node label should be the same
            df.iloc[i] = [int(node), node_attr['chr'], node_attr['chunk_start'], node_attr['chunk_end']]
            i += 1
        
        df.sort_values(by=['chr', 'chunk_start'], inplace=True) # sort according to chromosome, then chunk start
        df.reset_index(drop=True, inplace=True) # reset the sorted index
        df.to_csv(out_dir, index=False) # save to csv file
        time_elapsed = time.time()-start 
        print('time elapsed for exporting nodes: ', time_elapsed)


if __name__ == "__main__":
    args = get_args()
    data_dir = args.data_dir
    edge_dir = args.edge_dir
    node_dir = args.node_dir

    hic_graph = HicGraph(data_dir)
    hic_graph.report()
    hic_graph.export_graph(edge_dir, node_dir)

