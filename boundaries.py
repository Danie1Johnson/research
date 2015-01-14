"""
NAME_boundary(x) returns True if x is within the boundary and False otherwise.
"""
import numpy as np
import triangle_intersection as ti


def get_boundary(boundary_name, kwargs={}):
    try:
        b_fun = globals()[boundary_name+"_boundary"]
        boundary = lambda x: b_fun(x,**kwargs)
        return boundary  #, boundary_normal
    except KeyError, NameError:
        try:
            # For Building Game intermediates. Denoted 'polyname'.
            int_num = kwargs['int_num']
            n, dim, q0, masses, links, lengths, faces = bga.load_bg_int(manifold_name, int_num)
            boundary = lambda x: nonintersection_boundary(x, faces)
            return boundary  #, boundary_normal
       except ValueError, IndexError:
            pass
        print "ERROR:", manifold_name, "not found."
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
            if ti.triangle_intersection(t1.get_face(faces[j]), ti.get_face(faces[k])) = True:
                return False
    return True
