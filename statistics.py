import numpy as np

import bga_4_0 as bga

def get_statistic(stat_name, kwargs={}):
    try:
        s_ globals()[stat_name]
        return statistic
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


def test_1(x)

        
def angle_between_edges(x, v0, v1, v2):
    """
    Return the angle between edge (v0 v1) and edge (v2 v1). 
    """

    ans = np.dot(x[3*v0:3*v0+3] - x[3*v1:3*v1+3],
                 x[3*v2:3*v2+3] - x[3*v1:3*v1+3])
    ans /= numpy.linalg.norm(x[3*v0:3*v0+3] -x[3*v1:3*v1+3])
    ans /= numpy.linalg.norm(x[3*v2:3*v2+3] -x[3*v1:3*v1+3])
    return np.arccos(ans)

def dihedral_angle(x, v0, v1a, v1b, v2):
    """
    Return the angle between v0, the v1a--v1b midpoint and v2 
    """

    ans = np.dot(x[3*v0:3*v0+3] - 0.5*(x[3*v1a:3*v1a+3] + x[3*v1b:3*v1b+3]),
                 x[3*v2:3*v2+3] - 0.5*(x[3*v1a:3*v1a+3] + x[3*v1b:3*v1b+3]))
    ans /= numpy.linalg.norm(x[3*v0:3*v0+3] -0.5*(x[3*v1a:3*v1a+3] + x[3*v1b:3*v1b+3]))
    ans /= numpy.linalg.norm(x[3*v2:3*v2+3] -0.5*(x[3*v1a:3*v1a+3] + x[3*v1b:3*v1b+3])))

    return np.arccos(ans)
