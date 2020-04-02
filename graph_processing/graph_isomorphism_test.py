"""
Test if the loaded graph from h5 and csv is identical to the original gexf.
"""
import argparse
import networkx as nx
import networkx.algorithms.isomorphism as iso

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


if __name__=="__main__":
    args = get_args()
    data_old = args.data_old
    data_new = args.data_new
    graph_old = nx.read_gexf(data_old)
    graph_new = nx.read_gexf(data_new)
    
    print('matching topology...')
    print(nx.is_isomorphic(graph_old, graph_new))
    
    print('matching topology with node attributes...')
    nm = iso.categorical_node_match('chunk_start', None) # set default value to None
    print(nx.is_isomorphic(graph_old, graph_new, node_match=nm))
    nm = iso.categorical_node_match('chr', None) # set default value to None
    print(nx.is_isomorphic(graph_old, graph_new, node_match=nm))
    nm = iso.categorical_node_match('chunk_end', None) # set default value to None
    print(nx.is_isomorphic(graph_old, graph_new, node_match=nm))

    print('matching topology with edge attributes...')
    em = iso.numerical_edge_match('contactCount', None) # set default value to None
    print(nx.is_isomorphic(graph_old, graph_new, edge_match=em))
    em = iso.numerical_edge_match('p-value', None) # set default value to None
    print(nx.is_isomorphic(graph_old, graph_new, edge_match=em))
    em = iso.numerical_edge_match('q-value', None) # set default value to None
    print(nx.is_isomorphic(graph_old, graph_new, edge_match=em, node_match=nm))
    


    