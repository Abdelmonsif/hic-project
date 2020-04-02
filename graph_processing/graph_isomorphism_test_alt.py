"""
Alternative implementation of the graph isomorphism test by comparing the sorted node table and edge table.
"""
import argparse
import networkx as nx
import pandas as pd

def get_args():
    parser = argparse.ArgumentParser('python')

    parser.add_argument('-data_old',
                        default='../../test_data/test_4.gexf',
                        required=False,
                        help='directory of original hi-c graph')  

    parser.add_argument('-data_new',
                        default='../../test_data/test_4_recovered.gexf',
                        required=False,
                        help='directory of new hi-c graph')  
 
    return parser.parse_args()

def load_graph(graph_dir):
    """
    load the graph gexf file and return the sorted node table and edge table
    """
    g = nx.read_gexf(graph_dir)
    num_nodes = nx.number_of_nodes(g)
    num_edges = len(list(g.edges))
    
    '''form node table'''
    index = range(num_nodes)
    columns = ['node_id', 'chr', 'chunk_start', 'chunk_end']
    nodes = pd.DataFrame(index=index, columns=columns) # node table
    node_list = list(g.nodes)
    i = 0
    for node in node_list:
        node_attr = g.nodes[node]
        assert(node_attr['label']==node) # node id and node label should be the same
        nodes.iloc[i] = [int(node), int(node_attr['chr']), int(node_attr['chunk_start']), int(node_attr['chunk_end'])]
        i += 1
    nodes = nodes.sort_values(by=['chr', 'chunk_start'], inplace=False) # sort according to chromosome, then chunk start

    '''form edge table'''
    columns = ['source', 'target', 'contactCount', 'p_value', 'q_value']
    edge_table = pd.DataFrame(columns=columns)
    edge_list = list(g.edges)
    for edge in edge_list: # this loop can be parallelized in the future  
        edge_attr = g.edges[edge[0], edge[1]] # retrieve edge attribute   
        edge = tuple(map(int, edge)) # convert edge from string to integer for storage efficiency
        source = min(edge)
        target = max(edge)
        contactCount = edge_attr['contactCount']
        p_value = edge_attr['p-value']
        q_value = edge_attr['q-value']
        edge_row = pd.Series([source, target, contactCount, p_value, q_value], index=columns)
        edge_table = edge_table.append(edge_row, ignore_index=True)
    edge_table = edge_table.sort_values(by=['source', 'target'], inplace=False)
    nodes = nodes.reset_index(drop=True) # we don't compare indexes
    edge_table = edge_table.reset_index(drop=True) # we don't compare indexes
    print(edge_table)
    print(nodes)
    return nodes, edge_table


if __name__=="__main__":
    args = get_args()
    data_old = args.data_old
    data_new = args.data_new
    nodes_old, edge_table_old = load_graph(data_old)
    nodes_new, edge_table_new = load_graph(data_new)
    print('matching nodes......')
    if (nodes_old.equals(nodes_new)):
        print('passed!')
    else:
        print('failed!')
    print('matching edges......')
    if (edge_table_old.equals(edge_table_new)):
        print('passed!')
    else:
        print('failed!')