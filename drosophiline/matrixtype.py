"""All of the vectors and matrices are really numpy arrays.
"""

import numpy as np

def is_stochastic_vector(v):
    if len(v.shape) != 1:
        return False
    if any(x<0 for x in v):
        return False
    if abs(1.0 - v.sum()) > 1e-7:
        return False
    return True

def is_square_matrix(M):
    if len(M.shape) != 2:
        return False
    if len(set(M.shape)) != 1:
        return False
    return True

def is_right_stochastic_matrix(M):
    if not is_square_matrix(M):
        return False
    if not all(is_stochastic_vector(v) for v in M):
        return False
    return True
