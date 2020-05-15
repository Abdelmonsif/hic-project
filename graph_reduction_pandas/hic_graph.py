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
import itertools
from statistics import median

class HicGraph:
    def __init__(self, edge_dir, node_dir, snps_dir, write_gexf_dir, verbose):
        self.edge_dir = edge_dir
        self.node_dir = node_dir
        self.snps_dir = snps_dir
        self.write_gexf_dir = write_gexf_dir
        self.verbose = verbose
        

    def load_graph(self):
        """
        Load the original main graph. All edge attributes and each node's chromosome, chunk_start and chun_end are loaded.
        The graph reduction algorithm is going to work on the 2 loaded data structures: 'nodes' and 'edge_list'
        """
        '''Load the json file containing the SNP mapping.'''
        with open(self.snps_dir) as f:
            self.snp_map = json.load(f)

        self.edge_list, self.contactCount, self.p_values, self.q_values, self.edge_ids = self.__load_edge(self.edge_dir)
        
        self.nodes = self.__load_node(self.node_dir)
        self.node_list = self.nodes['node_id'] # get node list
        self.nodes.set_index('node_id', inplace=True)
        
        '''Generate a set of node ids. It is used to intersect with the set of node ids with SNPs.
        Might be deprecated in the future when using the whole main graph.'''
        self.node_id_set = set(self.nodes.index.values.tolist())
        
        self.edge_list = list(map(tuple, self.edge_list)) # edge list 

        self.edge_table = pd.DataFrame() # edge table       
        self.edge_table['source'] = [x[0] for x in self.edge_list]
        self.edge_table['target'] = [x[1] for x in self.edge_list]
        self.edge_table['contactCount'] = self.contactCount
        self.edge_table['p-value'] = self.p_values
        self.edge_table['q-value'] = self.q_values
        

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


    def load_patient(self, patient_dir):
        """
        Load the csv file containing SNPs of a patient, then add the locations of 
        SNPs to the nodes dataframe.
        """
        self.nodes['has_snp'] = False # add a column to indicate presence of SNPs
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
        print('computing nodes to merge...')
        start_time=time.time()
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
        to_merge = list(filter(None, to_merge)) # remove empty lists

        if self.verbose == 1:
            print('nodes to merge:')
            print(to_merge)

        self.__report_elapsed_time(start_time)

        self.__merge_nodes_1(to_merge)


    def __merge_nodes_1(self, to_merge):
        """
        Used by methods graph_reduce_1 and graph_reduce_2 to merge the specified nodes (non-SNP) into one node.
        """
        '''convert everything in edge table to list for future merging'''
        self.edge_table['source'] = [str([x]) for x in list(self.edge_table['source'])]
        self.edge_table['target'] = [str([x]) for x in list(self.edge_table['target'])]
        self.edge_table['contactCount'] = [str([x]) for x in list(self.edge_table['contactCount'])]
        self.edge_table['p-value'] = [str([x]) for x in list(self.edge_table['p-value'])]
        self.edge_table['q-value'] = [str([x]) for x in list(self.edge_table['q-value'])]

        if self.verbose == 1:
            print('nodes before reduction:')
            print(self.nodes)
            print('edge table before reduction')
            print(self.edge_table)
            print('merging nodes...')

        print('merging nodes...')
        start=time.time()
        '''merge nodes'''
        nodes_reduced = self.nodes[self.nodes['has_snp']==True].copy() # SNP-nodes, use .copy() to avoid SettingWithCopyWarning
        nodes_reduced['node_id'] = [str([idx]) for idx in list(nodes_reduced.index)] # convert single int to list
        nodes_reduced = nodes_reduced.set_index('node_id') # convert single int to list
        for node_list in to_merge:
            node_id = str(node_list)
            chr = self.nodes.loc[node_list[0], 'chr']
            chunk_start = self.nodes.loc[node_list[0], 'chunk_start']
            chunk_end = self.nodes.loc[node_list[-1], 'chunk_end']
            has_snp = False
            nodes_reduced.loc[node_id]= [chr, chunk_start, chunk_end, has_snp]
        nodes_reduced = nodes_reduced.sort_values(by=['chr', 'chunk_start'], inplace=False) # sort according to chromosome, then chunk start
        self.__report_elapsed_time(start)
        
        if self.verbose == 1:
            print('reduced nodes:')
            print(nodes_reduced)
            print('merging edges...')
        
        print('merging edges...')
        start=time.time()
        '''merge edges and edge attributes'''
        edges_reduced = self.edge_table.copy()
        for node_id, node in nodes_reduced.iterrows(): # traverse through the chromosomes
            if node['has_snp'] == False and len(eval(node_id)) > 1: # condidtion for merging
                node_id = eval(node_id) # convert from string to list
                old_edges = pd.DataFrame()
                for node in node_id:
                    source_edges = edges_reduced[edges_reduced['source']==str([node])].copy() # find the edges connected to this node (as source)
                    target_edges = edges_reduced[edges_reduced['target']==str([node])].copy() # find the edges connected to this node (as target)
                    target_edges.rename(columns={'source': 'target', 'target':'source'}, inplace=True) # exchange source and target 
                    target_edges = target_edges[['source', 'target', 'contactCount', 'p-value', 'q-value']] # move columns to make columns consistent        
                    old_edges = old_edges.append(source_edges)
                    old_edges = old_edges.append(target_edges)                        
                    index_to_remove = edges_reduced[(edges_reduced['source']==str([node])) | (edges_reduced['target']==str([node]))].index  # remove old edges
                    edges_reduced.drop(index_to_remove , inplace=True) # remove old edges
                
                if self.verbose == 1:
                    print('old_edges:\n', old_edges)
                target_nodes = list(set(list(old_edges['target']))) # list of targets (unique)
                target_nodes = [eval(x) for x in target_nodes] # list of targets (unique)
                if self.verbose == 1:
                    print('target_nodes: ', target_nodes)
                    print('node_id: ', node_id)
                    print('reduced nodes:', nodes_reduced)
                for target_node in target_nodes:
                    if not set(target_node).issubset(node_id): # avoid self-loops after merging
                        target_node_old_edges = old_edges[old_edges['target']==str(target_node)]
                        merged_edge = pd.Series([str(node_id), # source
                                                 str(target_node), # target
                                                 str(list((itertools.chain.from_iterable([eval(x) for x in list(target_node_old_edges['contactCount'])])))), 
                                                 str(list((itertools.chain.from_iterable([eval(x) for x in list(target_node_old_edges['p-value'])])))), 
                                                 str(list((itertools.chain.from_iterable([eval(x) for x in list(target_node_old_edges['q-value'])]))))], 
                                                 index=edges_reduced.columns) # new edge 
                        edges_reduced = edges_reduced.append(merged_edge, ignore_index=True)
                        if self.verbose == 1:                
                            print('merged edge:\n', merged_edge[['source', 'target', 'contactCount', 'p-value']])
                            print('edge table:\n', edges_reduced[['source', 'target', 'contactCount', 'p-value']])
                if self.verbose == 1:
                    print('----------------------------------------------------------------------')
        
        self.__report_elapsed_time(start)
        if self.verbose == 1:
            print('reduced nodes:\n', nodes_reduced)
            print('edge table:\n', edges_reduced)

        print('computing median...')
        start=time.time()
        '''compute median/mean after all merging done'''
        edges_reduced['contactCount'] = self.__compute_median(edges_reduced['contactCount']) # median of contactCount
        edges_reduced['p-value'] = self.__compute_median(edges_reduced['p-value']) # median of p-value
        edges_reduced['q-value'] = self.__compute_median(edges_reduced['q-value']) # median of q-value
        self.__report_elapsed_time(start)

        print('renaming nodes...')
        start=time.time()
        '''create node names with chr-chunk_start_chunk_end'''
        nodes_reduced = nodes_reduced.astype(str)
        chunk_start_series = nodes_reduced['chunk_start'].copy()
        chunk_end_series = nodes_reduced['chunk_end'].copy()
        nodes_reduced['node_names'] = nodes_reduced['chr'].str.cat(chunk_start_series, sep='-')
        nodes_reduced['node_names'] = nodes_reduced['node_names'].str.cat(chunk_end_series, sep='-')
        nodes_reduced.drop(columns=['chr', 'chunk_start', 'chunk_end'], inplace=True)
        node_name_dict = nodes_reduced['node_names'].to_dict() # dictionary to map original node names from node list to chr-chunk_start-chunk-end
        nodes_reduced.reset_index(inplace=True) # convert original node names from index to column
        nodes_reduced.rename(columns={"node_id": "original_nodes"}, inplace=True) # rename this column
        nodes_reduced.set_index('node_names', inplace=True)# set chr-chunk-start-chunk_end as new index
        self.nodes_reduced = nodes_reduced
        self.__report_elapsed_time(start)

        print('renaming edges...')
        start=time.time()
        '''rename source and target in edge table with chromosome-chunk_start-chunk_end'''
        edges_reduced['source'] = edges_reduced['source'].apply(lambda x: node_name_dict[x])
        edges_reduced['target'] = edges_reduced['target'].apply(lambda x: node_name_dict[x])
        self.edges_reduced = edges_reduced
        self.__report_elapsed_time(start)


    def __compute_median(self, edge_series):
        """
        Compute the median value for each edge attributes.
        """
        edge_series = [eval(x) for x in list(edge_series)]
        edge_series = [float(median(x)) for x in edge_series]
        return edge_series


    def graph_reduce_2(self):
        """
        A more efficient implementation of 'graph_reduce_1'
        """        
        '''compute nodes to merge'''
        print('computing nodes to merge...')
        start_time=time.time()
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
        to_merge = list(filter(lambda x: len(x)>1, to_merge)) # remove empty lists and lists with single elements
        if self.verbose == 1:
            print('nodes to merge:')
            print(to_merge)
        self.__report_elapsed_time(start_time)

        print('merging nodes...')
        start_time=time.time()
        '''creating a new column indicating node ids after merging'''
        self.nodes['new_id'] = self.nodes.index # need to execute this line everytime a new patient is loaded
        for node_list in to_merge: # compute new id by replacing old ones with the first id of nodes to merge
            self.nodes.loc[node_list, ['new_id']] = node_list[0]

        '''merge the nodes'''
        merged_nodes_set = list(set(list(self.nodes['new_id']))) # set of new node ids
        merged_nodes_list = [] # list of merged node series
        for node in merged_nodes_set:
            old_nodes = self.nodes[self.nodes['new_id']==node]
            merged_node_id = old_nodes.at[old_nodes.index[0],'new_id'] # chromosome
            chromosome = old_nodes.at[old_nodes.index[0],'chr'] # chromosome
            chunk_start = old_nodes['chunk_start'].min() # chunk start
            chunk_end = old_nodes['chunk_end'].max() # chunk end
            has_snp = old_nodes.at[old_nodes.index[0],'has_snp'] # has SNP?
            merged_node = pd.Series([merged_node_id, chromosome, chunk_start, chunk_end, has_snp], index=['node_id', 'chr','chunk_start','chunk_end','has_snp'])
            merged_nodes_list.append(merged_node)
        self.nodes_reduced = pd.concat(merged_nodes_list, axis=1).T # create new dataframe for final merged node table
        self.nodes_reduced.set_index('node_id', inplace=True) # set new index as the index of node table

        if self.verbose == 1:
            print(self.nodes_reduced)
        self.__report_elapsed_time(start_time)
        
        print('merging edges...')
        start_time=time.time()
        '''merge edges and generate a new edge table'''
        edge_table = self.edge_table.copy() # copy the original edge table, local to this function
        node_name_dict = self.nodes['new_id'].to_dict() # dictonry to map old nodes to new node ids
        edge_table['source'] = edge_table['source'].map(node_name_dict) # merging by renaming
        edge_table['target'] = edge_table['target'].map(node_name_dict) # merging by renaming
        source_target = edge_table[['source', 'target']].to_numpy() # make source always <= target
        source_target.sort(axis=1) # make source always <= target
        edge_table[['source', 'target']] = source_target # make source always <= target
        self_loop_idx = edge_table[edge_table['source']==edge_table['target']].index # remove self-loops (rows with same source and target)
        edge_table.drop(self_loop_idx , inplace=True) # remove self-loops (rows with same source and target)
        source_target_pairs = list(set(tuple(x) for x in edge_table[['source', 'target']].to_numpy())) # set of source-target pairs (source always <= target)
        edges_list = [] # list to store the new merged edges
        for source_target_pair in source_target_pairs: # compute median and append to new data frame
            source = source_target_pair[0]
            target = source_target_pair[1]
            edges = edge_table[(edge_table['source']==source) & (edge_table['target']==target)]
            contactCount = edges['contactCount'].median() # take median of edge attributes
            p_value = edges['p-value'].median() # take median of edge attributes
            q_value = edges['q-value'].median() # take median of edge attributes
            merged_edge = pd.Series([source, target, contactCount, p_value, q_value], index=edges.columns)
            edges_list.append(merged_edge)
        self.edges_reduced = pd.concat(edges_list, axis=1).T # reduced edge table
        convert_dict = {'source': int, 'target': int} 
        self.edges_reduced = self.edges_reduced.astype(convert_dict) 
        
        if self.verbose == 1:
            print('reduced edges:')
            print(self.edges_reduced)

        #print(self.edges_reduced.shape)
        self.__report_elapsed_time(start_time)

   
    def export_reduced_graph(self, reduced_node_dir, reduced_edge_dir):
        self.nodes_reduced.to_csv(reduced_node_dir, index=True) # save nodes
        self.edges_reduced.to_csv(reduced_edge_dir, index=False) # save edges


    def export_reduced_gexf(self, reduced_gexf_dir, reduced_graph_statistics_dir):
        """
        Export the reduced graph as gexf file.
        """
        reduced_graph = nx.Graph()

        '''nodes'''
        node_list = list(self.nodes_reduced.index.values) # node list
        reduced_graph.add_nodes_from(node_list) # add nodes to graph
        node_attr = self.nodes_reduced.to_dict('index') # make node dictionary 
        nx.set_node_attributes(reduced_graph, node_attr) # add node dictionary as node attributes
        
        '''edges'''
        self.edges_reduced['edge_names'] = list(zip(self.edges_reduced.source, self.edges_reduced.target)) # add reduced edge list to reduced edge table
        self.edges_reduced.drop(columns=['source', 'target'], inplace=True)# drop source and target column
        reduced_graph.add_edges_from(list(self.edges_reduced['edge_names'])) # add edges to graph
        self.edges_reduced.set_index('edge_names', inplace=True)# make new edge names index
        edge_attr = self.edges_reduced.to_dict('index') # make dictionary of edge attributes    
        nx.set_edge_attributes(reduced_graph, edge_attr)
        
        nx.write_gexf(reduced_graph, reduced_gexf_dir) # export the reduced graph as gexf file.

        '''statistics'''
        num_nodes_reduced = reduced_graph.number_of_nodes() # number of nodes
        num_edges_reduced = reduced_graph.number_of_edges() # number of edges
        avg_degree = sum([d[1] for d in reduced_graph.degree()])/num_nodes_reduced # average degree
        num_isolated = len(list(nx.isolates(reduced_graph))) # number of isolated nodes 
        density = nx.classes.function.density(reduced_graph) # density
        
        statistics_dict = {'num_nodes_reduced':num_nodes_reduced,
                           'num_edges_reduced':num_edges_reduced,
                           'avg_degree':avg_degree,
                           'num_isolated':num_isolated,
                           'density':density}

        with open(reduced_graph_statistics_dir, 'w') as outfile:
            json.dump(statistics_dict, outfile)

        #print('number of nodes:', num_nodes_reduced)
        #print('number of edges:', num_edges_reduced)
        #print('average degree:', avg_degree)
        #print('number of isolated nodes:', num_isolated)
        #print('density of reduced graph:', density)
        

        
        

    

