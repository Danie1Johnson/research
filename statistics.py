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
        boundary = lambda x: nonintersection_boundary(x, faces)
        return boundary  #, boundary_normal
    except ValueError, IndexError:
        pass
