"""
NAME_boundary(x) returns a numpy array with an entry for each boundary. 
Each entry negative for being within the boundary and positive for outside the boundary.

NAME_boundary_normal(x, b_num) Assuming x is on the b_num-th boundary (or very close to it),
return a unitvector (of the smae length as x) that is inwardly normal to the b_num-th boundary at x.
"""
import numpy as np

import bga_4_0 as bga

def get_manifold(manifold_name, kwargs={}):
    """
    Use manifold name to find and return constraint functions.
    """
    if manifold_name == None:
        return lambda x: np.empty(shape=(0)), np.empty(shape=(0,0)) 
    try:
        return globals()[manifold_name](**kwargs)
    except (KeyError, NameError):
        raise Exception("ERROR: " + manifold_name + " not found.")


###--------------------------------------------------------------------------
def building_game(poly_name=None, 
                  int_num=None, 
                  fixed_face=None,
                  fixed_com=False,
                  masses = None,
                  dim=3):
    """
    Load Building Game info for specified intermediate and return constraint functions.
    """
    try:
        q0, links, lengths, faces = bga.load_bg_int(poly_name, int_num)
    except:
        raise Exception("ERROR: Building game intermediate " + str(int_num) + 
                        " for " + polyname + " not found.")
    fixed_inds = [] 
    fixed_vals = []
                  
    if fixed_face != None:
        for j, vert in enumerate(faces[fixed_face]):
            # Only fix 6 dofs (i.e. dont over constrain)
            fixed_inds += [3*vert + k for k in range(3-j)]
            fixed_vals += [q0[3*vert+k] for k in range(3-j)]

    return linkage(links=links, 
                   lengths=lengths, 
                   fixed_inds=fixed_inds, 
                   fixed_vals=fixed_vals, 
                   fixed_com=fixed_com, 
                   masses=masses, 
                   dim=dim)
###--------------------------------------------------------------------------
def sphere(r=1.0):
    c = lambda x: sphere_c(x, r=r)
    C = lambda x: sphere_C(x, r=r)
    return c, C

def sphere_c(x, r=1.0):
    return ellipse_c(x, r=r)
            
def sphere_C(x):
    return ellipse_C(x)
###--------------------------------------------------------------------------
def ellipse(r=1.0, a=None, A=None):
    c = lambda x: ellipse_c(x, r=r, a=a, A=A)
    C = lambda x: ellipse_C(x, r=r, a=a, A=A)
    return c, C

def ellipse_c(x, r=1.0, a=None, A=None):
    c = []
    if a != None:
        c.append(np.dot(x/a,x) - r**2)
    elif A != None:
        c.append(np.dot(np.dot(x.T,A),x) - r**2)
    else:
        c.append(np.dot(x,x) - r**2)
    return np.array(c)

def ellipse_C(x, a=None, A=None):
    C = []
    if a != None:
        C.append(2.0*x/a)
    elif A != None:
        C.append(2.0*np.dot(A,x))
    else:
        C.append(2.0*x)
    return np.array(C)
###--------------------------------------------------------------------------            
def linkage(links=[], lengths=[], fixed_inds=[], fixed_vals=[], fixed_com=False, masses=None, dim=3):
    c = lambda x: linkage_c(x, 
                            links=links, 
                            lengths=lengths, 
                            fixed_inds=fixed_inds, 
                            fixed_vals=fixed_vals, 
                            fixed_com=fixed_com, 
                            masses=masses, 
                            dim=dim)
    C = lambda x: linkage_C(x, 
                            links=links,
                            fixed_inds=fixed_inds, 
                            fixed_com=fixed_com, 
                            masses=masses, 
                            dim=dim)
    return c, C

def linkage_c(q, links=[], lengths=[], fixed_inds=[], fixed_vals=[], fixed_com=False, masses=None, dim=3):
    """
    Return np array of length m containting each constraint evaluated at q.
    """
    if len(fixed_inds) != len(fixed_vals):
        raise Exception("ERROR: fixed_inds and fixed_vals lengths must correspond")
    nf = len(fixed_inds)
    if fixed_com == True:
        ncom = dim
        if masses == None:
            masses = np.ones(len(q)/dim)
        if nf != 0:
            raise Exception("Error: cannot both fix face and center of mass")
    else:
        ncom = 0
    m = len(links) + nf + ncom
    c = np.zeros((m,))

    # Fixed faces
    for k in range(nf):
        c[k] = q[fixed_inds[k]] - fixed_vals[k]

    # Fixed center of mass
    for j in range(ncom):
        c[j] = np.dot(q[j::dim], masses)

    # Link constraints
    for i, link in enumerate(links):
        #print i, link, len(link), len(lengths)
        c[i+nf+ncom] = sum((q[dim*link[0]:dim*link[0] + dim] - 
                            q[dim*link[1]:dim*link[1] + dim])**2) - lengths[i]**2
    return c

def linkage_C(q, links, fixed_inds=[], fixed_com=False, masses=None, dim=3):
    """
    Compute Jacobian matrix of c at q.
    """    
    nf = len(fixed_inds)
    if fixed_com == True:
        ncom = dim
        if masses == None:
            masses = np.ones(len(q)/dim)
        if nf != 0:
            raise Exception("Error: cannot both fix face and center of mass")
    else:
        ncom = 0
    m = len(links) + nf + ncom
    n = len(q)                                    
    C = np.zeros((m,n))
    # Fixed faces
    for i in range(nf):
        C[i,fixed_inds[i]] = 1.0

    # Fixed center of mass
    for j in range(ncom):
        C[j,j::dim] = masses

    # Link constraints
    for k, link in enumerate(links):
        for d in range(dim):
            C[k+nf+ncom,link[0]*dim + d] += 2.0*(q[link[0]*dim + d] - q[link[1]*dim + d])
            C[k+nf+ncom,link[1]*dim + d] += -2.0*(q[link[0]*dim + d] - q[link[1]*dim + d])
    return C
