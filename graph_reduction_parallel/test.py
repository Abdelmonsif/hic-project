"""
Test how to delete rows from a numpy array efficiently
"""
import numpy as np


if __name__ == "__main__":
    a = np.array(list(range(100)))
    a = np.reshape(a, (10,10)) 
    print(a)

    to_remove = [20, 30, 50, 90]
    l = []
    for x in to_remove:
         idx = np.nonzero(a[:,0]==x)
         l.append(idx[0].item())
    l = np.array(l)
    print(l)
    a = np.delete(a, l, 0)
    print(a)
