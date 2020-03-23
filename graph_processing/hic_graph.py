"""
The hi-c graph class containing all methods needed for graph reduction
"""
import h5py
import numpy as np
import networkx as nx
import pandas as pd
from os import listdir
from os.path import isfile, join
import time
import json

class HicGraph:
    def __init__(self, edge_dir, node_dir, snps_dir, write_gexf_dir):
        self.edge_dir = edge_dir
        self.node_dir = node_dir
        self.snps_dir = snps_dir
        self.write_gexf_dir = write_gexf_dir
        

    def load_graph(self):
        """
        Load the original main graph. All edge attributes and each node's chromosome, chunk_start and chun_end are loaded.
        The graph reduction algorithm is going to work on the 2 loaded data structures: 'nodes' and 'edge_list'
        """
        print('loading the main graph...')
        start=time.time()

        self.__load_snp_map() # load SNP map

        self.edge_list, self.contactCount, self.p_values, self.q_values, self.edge_ids = self.__load_edge(self.edge_dir)
        
        self.nodes = self.__load_node(self.node_dir)
        self.node_list = self.nodes['node_id'] # get node list
        self.nodes.set_index('node_id', inplace=True)
        self.nodes['has_snp'] = False # add a column to indicate presence of SNPs
        
        # Generate a set of node ids. It is used to intersect with the set of node ids with SNPs.
        # Might be deprecated in the future when using the whole main graph.
        self.node_id_set = set(self.nodes.index.values.tolist())
        
        self.edge_list = list(map(tuple, self.edge_list))
        
        print('loading finished. Time for loading the graph:')
        self.__report_elapsed_time(start)


    def __report_elapsed_time(self, start):
        end = time.time()   
        time_elapsed = end - start
        hour = time_elapsed//3600
        time_elapsed = time_elapsed - hour * 3600
        minute = time_elapsed//60
        time_elapsed = time_elapsed - minute * 60
        print('{}:{}:{}'.format(int(hour), int(minute), round(time_elapsed)))


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
        Load nodes from the csv file produced by graph_preprocess.py
        """
        nodes = pd.read_csv(node_dir)
        return nodes


    def __load_snp_map(self):
        """
        Load the json file containing the SNP mapping.
        """
        with open(self.snps_dir) as f:
            self.snp_map = json.load(f)


    def load_patient(self, patient_dir):
        """
        Load the csv file containing SNPs of a patient, then add the locations of 
        SNPs to the nodes dataframe.
        """
        patient_snp = pd.read_csv(patient_dir, sep='	') # load patient SNPs as a dataframe
        self.snp_cols = [] # list containing all the SNPs of the patient
        snp_cols_1 = patient_snp.columns[(patient_snp == 1).iloc[0]].tolist()
        snp_cols_2 = patient_snp.columns[(patient_snp == 2).iloc[0]].tolist()
        self.snp_cols.extend(snp_cols_1)
        self.snp_cols.extend(snp_cols_2)
        snp_locations = [] # find the locations (node ids) of the snps
        num_missing_snp = 0 # number of missing snp for this patient, ignore them
        for snp in self.snp_cols:
            try:
                snp_locations.append(self.snp_map[snp][-1]) # last element is node id
                #print(self.snp_map[snp])
            except:
                num_missing_snp += 1
        snp_locations = set(snp_locations) # there are multiple SNPs on a single node, so take the set of this list to remove duplicates
        snp_locations = list(snp_locations.intersection(self.node_id_set)) # this step may be redundant when using main graph as input
        self.nodes.loc[snp_locations, 'has_snp'] = True # True if there is SNP, False if no SNP


    def export_original_graph(self):
        """
        Write the loaded graph to gexf file. Need to call load_graph(self) first.
        The goal of this method is to verify the correctness of loading by graph isomorphism test.
        """
        hic_graph = nx.Graph() # create empty graph
        hic_graph.add_nodes_from(self.node_list.tolist()) # add nodes
        node_attr = self.nodes.to_dict('index') # make node dictionary 
        nx.set_node_attributes(hic_graph, node_attr) # add node dictionary as node attributes
        
        edge_attr = pd.DataFrame() # dataframe for edges
        edge_attr['contactCount'] = self.contactCount
        edge_attr['p-value'] = self.p_values
        edge_attr['q-value'] = self.q_values
        edge_attr['id'] = self.edge_ids 
        edge_attr['edge_list'] = self.edge_list
        edge_attr.set_index('edge_list', inplace=True) # set edge list as index
        edge_attr = edge_attr.to_dict('index') # convert to dictionary 
        hic_graph.add_edges_from(self.edge_list) # add edges
        nx.set_edge_attributes(hic_graph, edge_attr)
        
        nx.write_gexf(hic_graph, self.write_gexf_dir) # for testing graph isomorphism only
    

    def report(self):
        """
        Report info about this graph
        """
        #print('number of nodes in the graph:', nx.number_of_nodes(self.hic_graph))
        #print('attribute names of each node:', self.node_attribute_list)
        #print('attribute names of each edge:', self.edge_attribute_list)
        return None

    
    def graph_reduce_1(self):
        """
        Merge the nodes without SNP to nodes with SNP. The result graph will be a new NetworkX graph object.
        All non-SNP nodes are merged together.
        """
        to_merge = []        
        for chr in range(1,24): # for each chromosome
            rows_to_merge = [] # list of lists of rows in 23 chr dataframes to merge
            df_chr = self.nodes[self.nodes['chr'] == chr]
            snp_rows = np.where(df_chr.has_snp==True)[0] # row number of SNP nodes
            num_rows = df_chr.shape[0]
            rows = list(range(num_rows))

            if len(snp_rows) > 0:
                start = 0
                for snp_node in snp_rows: # generate lists of merged nodes
                    end = snp_node
                    if end > start: # avoid empty list when first node is SNP-node or consecutive SNP-nodes
                        rows_to_merge.append(rows[start:end]) # get nodes whose row_numbers are smaller than snp_node 
                    start = snp_node + 1 # next segment
                end = num_rows # last segment
                if end > start: # only do so when last node is not a SNP-node
                    rows_to_merge.append(rows[start:end]) # get nodes whose row_numbers are smaller than snp_node 
            elif len(snp_rows) == 0: # all nodes of this chromosome are non-SNP
                rows_to_merge.append(rows)
            
            node_ids = list(df_chr.index)
            node_ids_to_merge = [list(map(node_ids.__getitem__, rows))  for rows in rows_to_merge] # get the node ids of the nodes to merge into one node on this chromosome
            to_merge.extend(node_ids_to_merge)
            #print(df_chr)
            #print(rows_to_merge) 
            #print(node_ids_to_merge)
        to_merge = list(filter(None, to_merge)) # remove empty lists
        print('nodes to merge:')
        print(to_merge)
        
        self.__merge_nodes(to_merge)


    def graph_reduce_2(self):
        """
        Merge the non-SNP nodes to nearest SNP nodes according to the distances.
        """
        return None


    def __merge_nodes(self, to_merge):
        """
        Used by methods graph_reduce_1 and graph_reduce_2 to merge the specified nodes (non-SNP) into one node.
        """
        print(self.nodes)
        print(self.edge_list)
        reduced_graph = nx.Graph()
        # nodes
        # neighbors
        # edges
        # edge attributes
        

    

