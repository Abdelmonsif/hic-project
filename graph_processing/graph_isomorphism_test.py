"""
Test if the loaded graph from h5 and csv is identical to the original gexf.
"""
import argparse
import networkx as nx
import networkx.algorithms.isomorphism as iso

def get_args():
    parser = argparse.ArgumentParser('python')

    parser.add_argument('-data_old',
                        default='../../test_data/test_0.gexf',
                        required=False,
                        help='directory of original hi-c graph')  

    parser.add_argument('-data_new',
                        default='../../test_data/test_0_recovered.gexf',
                        required=False,
                        help='directory of new hi-c graph')  
 
    return parser.parse_args()


if __name__=="__main__":
    args = get_args()
    data_old = args.data_old
    data_new = args.data_new
    graph_old = nx.read_gexf(data_old)
    graph_new = nx.read_gexf(data_new)
    
    '''Test isomorphism with all combinations of node and edge attributes. 
       This way seems faster than matching nodes or edges along.'''
    node_attr_names = ['chr', 'chunk_start', 'chunk_end']
    edge_attr_names = ['contactCount', 'p-value', 'q-value']
    for node_attr in node_attr_names:
        for edge_attr in edge_attr_names:
            print('matching graphs with {} and {}...'.format(node_attr, edge_attr))
            nm = iso.categorical_node_match(node_attr, None) # set default value to None
            em = iso.numerical_edge_match(edge_attr, None) # set default value to None
            if nx.is_isomorphic(graph_old, graph_new, edge_match=em, node_match=nm):
                print('Pass!')
            else:
                print('Fail!')
            print('-------------------------------------------')



    