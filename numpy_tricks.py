import numpy as np

def expand_base_mones(items, size_inc):
    new_items = np.concatenate((items, -np.ones(size_inc, dtype=items.dtype)))
    return new_items

def expand_base_empty(items, size_inc):
    new_items = np.concatenate((items, np.empty(size_inc, dtype=items.dtype)))
    return new_items

def expand_base_zeros(items, size_inc):
    new_items = np.concatenate((items, np.zeros(size_inc, dtype=items.dtype)))
    return new_items


def expand_rows(array, row_inc):
    num_rows, num_cols = array.shape
    new_rows = np.empty((num_rows + row_inc, num_cols), dtype=array.dtype)
    new_rows[:num_rows, :] = array
    return new_rows

