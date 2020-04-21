"""
Parallel implementation of graph_reduction_np.py
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

