"""
Reduce the size of the graph by contracting edges.
"""
import argparse
import networkx as nx
from ast import literal_eval


def get_args():
    parser = argparse.ArgumentParser('python')

    parser.add_argument('-data_dir',
                        default='../../data/G_snps_23.gexf',
                        required=False,
                        help='directory of gene expression graph')  
                 
    return parser.parse_args()


class HicGraph:
    def __init__(self, data_dir):
        self.hic_graph = nx.read_gexf(data_dir) 
        self.node_attribute_list = ["chr",         
                                    "chunk_start", 
                                    "chunk_end",       
                                    "kr_norm",         
                                    "numberoftf",      
                                    "tfs",             
                                    "numberofmethy",   
                                    "methys"          
                                    "numberoflncrna",  
                                    "lncrna",          
                                    "numberofsnps",
                                    "snps",
                                    "numberofenhancers",
                                    "enhancers"]
        self.edge_attribute_list = ["contactCount", "p-value", "q-value"]    
        print(nx.degree(self.hic_graph))

    def report(self):
        """
        Report info about this graph
        """
        print('number of nodes in the graph:', nx.number_of_nodes(self.hic_graph))
        print('attribute names of each node:', self.node_attribute_list)
        print('attribute names of each edge:', self.edge_attribute_list)


    def node_list_gen(self):
        """
        Returns an iterator over the graph nodes
        """
        node_list = nx.nodes(self.hic_graph)
        node_list_int = [int(x) for x in node_list]
        node_list_int.sort()
        node_list = [str(x) for x in node_list_int]
        return node_list

    def get_node_attribute(self, node, attr_name):
        """
        node: integer, id of a node.
        attr_name: string, name of an attribute.
        """
        return self.hic_graph.nodes[node][attr_name]


    def neighbor_on_the_chain(self, node):
        """
        Returns a neighbor of the input node
        """

    def reduce_node(self, node):
        """
        Merge the input node with its neighbor.
        """

if __name__ == "__main__":
    args = get_args()
    data_dir = args.data_dir
    hic_graph = HicGraph(data_dir)
    hic_graph.report()
    node_list = hic_graph.node_list_gen()
    print(node_list)

    print('---------------------------------------------------------')
    test_attr = hic_graph.get_node_attribute(node_list[0], 'chr')
    print(test_attr)
    print(type(test_attr))

    print('---------------------------------------------------------')
    test_attr = hic_graph.get_node_attribute(node_list[0], 'chunk_start')
    print(test_attr)
    print(type(test_attr))

    print('---------------------------------------------------------')
    test_attr = hic_graph.get_node_attribute(node_list[0], 'chunk_end')
    print(test_attr)
    print(type(test_attr))

    print('---------------------------------------------------------')
    test_attr = hic_graph.get_node_attribute(node_list[0], 'kr_norm')
    print(test_attr)
    print(type(test_attr))
    test_attr = literal_eval(test_attr)
    print(test_attr)
    print(type(test_attr))

    print('---------------------------------------------------------')
    test_attr = hic_graph.get_node_attribute(node_list[0], 'numberoftf')
    print(test_attr)
    print(type(test_attr))

    print('---------------------------------------------------------')
    test_attr = hic_graph.get_node_attribute(node_list[0], 'tfs')
    print(test_attr)
    print(type(test_attr))
    test_attr = literal_eval(test_attr)
    print(test_attr)
    print(type(test_attr))

    print('---------------------------------------------------------')
    test_attr = hic_graph.get_node_attribute(node_list[0], 'numberofmethy')
    print(test_attr)
    print(type(test_attr))

    print('---------------------------------------------------------')
    test_attr = hic_graph.get_node_attribute(node_list[0], 'methys')
    print(test_attr)
    print(type(test_attr))
    test_attr = literal_eval(test_attr)
    print(test_attr)
    print(type(test_attr))