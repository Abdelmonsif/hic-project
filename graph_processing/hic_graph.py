"""
The hi-c graph class containing all methods needed for graph reduction
"""
import h5py
import numpy as np
import networkx as nx
import pandas as pd

class HicGraph:
    def __init__(self, edge_dir, node_dir):
        self.edge_dir = edge_dir
        self.node_dir = node_dir

    def load_graph(self):
        """
        Load the graph. All edge attributes and each node's chromosome, chunk_start and chun_end are loaded.
        """
        edge_list, contactCount, p_values, q_values, edge_ids = self.__load_edge(self.edge_dir)
        nodes = self.__load_node(self.node_dir)
        
        self.hic_graph = nx.Graph() # create empty graph
        node_list = nodes['node_id'] # get node list

        self.hic_graph.add_nodes_from(node_list.tolist()) # add nodes
        nodes.set_index('node_id', inplace=True)
        node_attr = nodes.to_dict('index') # make node dictionary 
        nx.set_node_attributes(self.hic_graph, node_attr) # add node dictionary as node attributes
        
        edge_list = list(map(tuple, edge_list))
        self.hic_graph.add_edges_from(edge_list) # add edges
        edge_attr = pd.DataFrame() # dataframe for edges
        edge_attr['contactCount'] = contactCount
        edge_attr['p-value'] = p_values
        edge_attr['q-value'] = q_values
        edge_attr['id'] = edge_ids 
        edge_attr['edge_list'] = edge_list
        edge_attr.set_index('edge_list', inplace=True) # set edge list as index
        edge_attr = edge_attr.to_dict('index') # convert to dictionary 
        nx.set_edge_attributes(self.hic_graph, edge_attr)

        nx.write_gexf(self.hic_graph, "./test_new.gexf")

    def __load_edge(self, edge_dir):
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

    def __load_node(self, node_dir):
        """
        Load nodes from the 
        """
        nodes = pd.read_csv(node_dir)
        return nodes

    def report(self):
        """
        Report info about this graph
        """
        print('number of nodes in the graph:', nx.number_of_nodes(self.hic_graph))
        print('attribute names of each node:', self.node_attribute_list)
        print('attribute names of each edge:', self.edge_attribute_list)