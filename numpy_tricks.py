import numpy as np

def expand(items, size_inc):
    new_items = np.concatenate((items, np.empty(size_inc, dtype=items.dtype)))
    return new_items