import networkx as nx

g = nx.Graph()
#nodes = [[1,2],3, 4,[5,6,7]]
nodes = [1,2,3, '[5,6,7]']
g.add_nodes_from(nodes) # add nodes
print(list(g.nodes))