"""
Return function boundary_name(x) that returns True if x is within the boundary and False otherwise.
"""
import numpy as np

import triangle_intersection as ti
import bga_4_0 as bga
import statistics as sts

def get_boundary(boundary_name, kwargs={}, binary=False):
    if boundary_name == None:
        if binary == True:
            return lambda x, y: True
        else:
            return lambda x: True
    return globals()[boundary_name](**kwargs)
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
def dihedrals(poly_name=None, int_num=None):
    """
    For building game intermediates, test for dihedral angle sign switches.
    """
    try:
        q0, links, lengths, faces = bga.load_bg_int(poly_name, int_num)
    except (ValueError, IndexError):
        raise Exception("ERROR: Building game intermediate " + str(int_num) + 
                " for " + polyname + " not found.")
    
    dihedral_inds = []
    F = len(faces)
    for j in range(F):
        for k in range(j+1, F):
            if adjacent_faces(j, k, faces) == True:
                dihedral_inds.append(order_verts(faces[j], faces[k]))
    
    return lambda x, y: dihedrals_fun(x, y, dihedral_inds=dihedral_inds)


def dihedrals_fun(x, y, dihedral_inds=None):
    """
    Compare the dihedrals of x and y and check for sign change. 
    """
    x_dihedrals = get_dihedrals(x, dihedral_inds)
    y_dihedrals = get_dihedrals(y, dihedral_inds)

    if max(abs(x_dihedrals - y_dihedrals)) > 2.0:
        return False
    else:
        return True

def get_dihedrals(y, dihedral_inds=None):
    """
    Compure the dihedral angles at y.
    """
    dihedrals = np.zeros(len(dihedral_inds))
    for k, di in enumerate(dihedral_inds):
        dihedrals[k] = sts.signed_dihedral_angle(y, di[0], di[1], di[2], di[3])
    return dihedrals

def order_verts(v1, v2):
    """
    get order of vert inds for dihedral calculation.
    """
    common_verts = set(v1).intersection(set(v2))
    ordered_verts = [list(set(v1).difference(common_verts))[0],
                     list(common_verts)[0],
                     list(common_verts)[1],
                     list(set(v2).difference(common_verts))[0]]
    return ordered_verts
###--------------------------------------------------------------------------
