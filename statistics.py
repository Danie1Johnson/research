import numpy as np
import numpy.linalg

import bga_4_0 as bga

def get_stat(stat_name, kwargs={}):
    try:
        s_fun = globals()[stat_name](**kwargs)
        return s_fun
    except KeyError, NameError:
        print "ERROR:", stat_name, "not found."
        raise


def bg_attachment(x, poly_name=None, int_num=None): 
    try:
        # For Building Game intermediates. Denoted 'polyname'.
        int_num = kwargs['int_num']
        q0, links, lengths, faces = bga.load_bg_int(boundary_name, int_num)
        V, E, F, S, species, f_types, adj_list, dual = bga.get_poly(poly_name)
        ints, ids, paths, shell_int, shell_paths, edges, shell_edge = bga.get_bg_ss(poly_name)
        try:
            int_faces = ints[int_num]
        except IndexError:
            print "ERROR:", poly_name, "does not have an intermediate", int_num
            raise
    except ValueError, IndexError:
        raise


def test_1():
    """
    Test for uniformity of dihedral angle in two linked triangles.
    """
    #return lambda x: np.array([dihedral_angle(x, 0, 2, 3, 1)]) 
    return lambda x: np.array([signed_dihedral_angle(x, 0, 2, 3, 1)]) 

def test_2():
    """
    Test for uniformity of dihedral angle in two linked triangles.
    """
    return lambda x: np.array([signed_dihedral_angle(x, 0, 3, 4, 1), 
                               signed_dihedral_angle(x, 3, 0, 4, 2)]) 
        
def angle_between_edges(x, v0, v1, v2):
    """
    Return the angle between edge (v0 v1) and edge (v2 v1). 
    """

    ans = np.dot(x[3*v0:3*v0+3] - x[3*v1:3*v1+3],
                 x[3*v2:3*v2+3] - x[3*v1:3*v1+3])
    ans /= numpy.linalg.norm(x[3*v0:3*v0+3] -x[3*v1:3*v1+3])
    ans /= numpy.linalg.norm(x[3*v2:3*v2+3] -x[3*v1:3*v1+3])
    return np.arccos(ans)

#def dihedral_angle(x, v0, v1a, v1b, v2):
#    """
#    Return the angle between v0, the v1a--v1b midpoint and v2 
#    """
#
#    #ans = np.dot(x[3*v0:3*v0+3] - 0.5*(x[3*v1a:3*v1a+3] + x[3*v1b:3*v1b+3]),
#    #             x[3*v2:3*v2+3] - 0.5*(x[3*v1a:3*v1a+3] + x[3*v1b:3*v1b+3]))
#    #ans /= numpy.linalg.norm(x[3*v0:3*v0+3] - 0.5*(x[3*v1a:3*v1a+3] + x[3*v1b:3*v1b+3]))
#    #ans /= numpy.linalg.norm(x[3*v2:3*v2+3] - 0.5*(x[3*v1a:3*v1a+3] + x[3*v1b:3*v1b+3]))
#    #return np.arccos(ans)


#def find_dihedrals(f1, f2, x, faces):
#    """
#    Take the index of two faces of a BG intermediate and compute the dihedral angle between them
#    """
#    n1 = triangle_normal(f1, x, faces)
#    n2 = triangle_normal(f2, x, faces)
#    return np.arccos(np.dot(n1, n2))
#
#def triangle_normal(f, x, faces):
#    """
#    Find normal vector to triangle. 
#    """
#    v1 = x[3*faces[f][0]:3*faces[f][0]+3] - x[3*faces[f][1]:3*faces[f][1]+3] 
#    v2 = x[3*faces[f][2]:3*faces[f][2]+3] - x[3*faces[f][1]:3*faces[f][1]+3] 
#    n = np.cross(v1, v2)
#    return n/numpy.linalg.norm(n)


def triangle_normal(x, va, vb, vc):
    """
    Find the normal vector at va
    """ 
    return np.cross(x[3*vb:3*vb+3] - x[3*va:3*va+3], x[3*vc:3*vc+3] - x[3*va:3*va+3]) 

def signed_dihedral_angle2(x, v0, v1a, v1b, v2):
    """
    Return the angle between v0, the v1a--v1b midpoint and v2 
    """
    # Triangle normals at v1a
    n1 = triangle_normal(x, v1a, v1b, v0)
    n2 = triangle_normal(x, v1a, v1b, v2)

    ans = np.arccos(np.dot(n1, n2)/(numpy.linalg.norm(n1)*numpy.linalg.norm(n2)))
    
    ###
    ans = np.pi - ans

    if np.dot(n1, x[3*v2:3*v2+3]- x[3*v1a:3*v1a+3]) > 0:
        return ans
    else:
        return 2.0*np.pi - ans

def signed_dihedral_angle(x, v0, v1a, v1b, v2):
    """
    Return the angle between v0, the v1a--v1b midpoint and v2 
    """

    ans = np.dot(x[3*v0:3*v0+3] - 0.5*(x[3*v1a:3*v1a+3] + x[3*v1b:3*v1b+3]),
                 x[3*v2:3*v2+3] - 0.5*(x[3*v1a:3*v1a+3] + x[3*v1b:3*v1b+3]))
    ans /= numpy.linalg.norm(x[3*v0:3*v0+3] -0.5*(x[3*v1a:3*v1a+3] + x[3*v1b:3*v1b+3]))
    ans /= numpy.linalg.norm(x[3*v2:3*v2+3] -0.5*(x[3*v1a:3*v1a+3] + x[3*v1b:3*v1b+3]))

    n = np.cross(x[3*v1a:3*v1a+3] - x[3*v0:3*v0+3], x[3*v1b:3*v1b+3] - x[3*v0:3*v0+3]) 
    theta = np.arccos(ans)
    if np.dot(x[3*v2:3*v2+3] - x[3*v0:3*v0+3], n) > 0.0:
        return theta
    else:
        return 2.0*np.pi - theta
    
