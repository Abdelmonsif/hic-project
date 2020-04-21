"""
Parallel implementation of the methods in hic_graph. 
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