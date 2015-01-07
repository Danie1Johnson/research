"""
NAME_boundary(x) returns a numpy array with an entry for each boundary. 
Each entry negative for being within the boundary and positive for outside the boundary.

NAME_boundary_normal(x, b_num) Assuming x is on the b_num-th boundary (or very close to it),
return a unitvector (of the smae length as x) that is inwardly normal to the b_num-th boundary at x.
"""
import numpy as np

def get_boundary(boundary_name, kwargs={}):
    try:
        b_fun = globals()[boundary_name+"_boundary"]
        bn_fun = globals()[boundary_name+"_boundary_normal"]
        boundary = lambda x: b_fun(x,**kwargs)
        boundary_normal = lambda x, y: bn_fun(x, y, **kwargs)
        return boundary, boundary_normal
    except KeyError, NameError:
        print "ERROR:", boundary_name, "not found."

###--------------------------------------------------------------------------
def positive_boundary(x):
    return -x

def positive_boundary_normal(x, b_num):
    n = np.zeros_like(x)
    n[b_num] = 1.0
    return n
###--------------------------------------------------------------------------
def none_boundary(x):
    return []

def none_boundary_normal(x, b_num):
    return None
###--------------------------------------------------------------------------
