"""
Return function boundary_name(x) that returns True if x is within the boundary and False otherwise.
"""
import triangle_intersection as ti
import bga_4_0 as bga

def get_boundary(boundary_name, kwargs={}):
    if boundary_name == None:
        return lambda x: True
    try:
        return globals()[boundary_name](**kwargs)
    except (KeyError, NameError):
        raise Exception("ERROR: " + boundary_name + " not found.")
        
###--------------------------------------------------------------------------
def positive_boundary():
    return lambda x: min(x) >= 0.0

####--------------------------------------------------------------------------
#def none():
#    return lambda x: True

###--------------------------------------------------------------------------
def self_intersection(poly_name=None, int_num=None):
    """
    For building game intermediates, test for intersection of triangles.
    """
    try:
        q0, links, lengths, faces = bga.load_bg_int(poly_name, int_num)
        return lambda x: self_intersection_fun(x, faces)
    except (ValueError, IndexError):
        raise Exception("ERROR: Building game intermediate " + str(int_num) + 
                " for " + polyname + " not found.")

def adjacent_faces(j, k, faces):
    return len(set(faces[j]).intersection(set(faces[k]))) == 2
 
def self_intersection_fun(x, faces):
    F = len(faces)
    for j in range(F):
        for k in range(j+1, F):
            if adjacent_faces(j, k, faces) == True:
                 continue
            if ti.triangle_intersection(ti.get_face(x, faces[j]), 
                                        ti.get_face(x, faces[k]), 
                                        scaling=0.99) == True:
                return False
    return True
###--------------------------------------------------------------------------
