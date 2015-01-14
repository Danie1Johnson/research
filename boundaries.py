"""
NAME_boundary(x) returns True if x is within the boundary and False otherwise.
"""
import numpy as np
import triangle_intersection as ti
import bga_4_0 as bga

def get_boundary(boundary_name, kwargs={}):
    try:
        b_fun = globals()[boundary_name+"_boundary"]
        boundary = lambda x: b_fun(x,**kwargs)
        return boundary  #, boundary_normal
    except KeyError, NameError:
        try:
            # For Building Game intermediates. Denoted 'polyname'.
            int_num = kwargs['int_num']
            n, dim, q0, masses, links, lengths, faces = bga.load_bg_int(boundary_name, int_num)
            boundary = lambda x: nonintersection_boundary(x, faces)
            return boundary  #, boundary_normal
        except ValueError, IndexError:
            pass
        print "ERROR:", boundary_name, "not found."
        raise

###--------------------------------------------------------------------------
def positive_boundary(x):
    return min(x) < 0.0

###--------------------------------------------------------------------------
def none_boundary(x):
    return True

###--------------------------------------------------------------------------
def nonintersection_boundary(x, faces):
    F = len(faces)
    for j in range(F):
        for k in range(j+1, F):
            if ti.triangle_intersection(ti.get_face(x, faces[j]), ti.get_face(x, faces[k])) == True:
                return False
    return True
