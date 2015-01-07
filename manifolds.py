"""
NAME_boundary(x) returns a numpy array with an entry for each boundary. 
Each entry negative for being within the boundary and positive for outside the boundary.

NAME_boundary_normal(x, b_num) Assuming x is on the b_num-th boundary (or very close to it),
return a unitvector (of the smae length as x) that is inwardly normal to the b_num-th boundary at x.
"""
import numpy as np

def get_manifold(manifold_name, kwargs={}):
    """
    Use manifold name to find and return constraint functions.
    """
    try:
        c_fun = globals()[manifold_name+"_c"]
        C_fun = globals()[manifold_name+"_C"]
        c = lambda x: c_fun(x,**kwargs)
        C = lambda x: C_fun(x,**kwargs)
        return c, C
    except KeyError, NameError:
        try:
            [poly_name, int_num] = manifold_name.split('__')
            n, dim, q0, masses, links, lengths, faces = bga.load_bg_int(poly_name, int_num)
            c = lambda x: linkage_c_fun(x, links, lengths)
            C = lambda x: linkage_C_fun(x, links)
            return c, C
        except ValueError, IndexError:
            pass
        print "ERROR:", manifold_name, "not found."
        raise

###--------------------------------------------------------------------------
def sphere_c(x, r=1.0):
    return ellipse_c(x, r=r)
            
def sphere_C(x):
    return ellipse_C(x)
###--------------------------------------------------------------------------
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
            
def linkage_c_fun(q, links, lengths, dim=3):
    """
    Return np array of length m containting each constraint evaluated at q.
    """
    m = len(links)
    c = np.zeros((m,))
    for i, link in enumerate(links):
        c[i] = sum((q[dim*link[0]:dim*link[0] + dim] - 
                    q[dim*link[1]:dim*link[1] + dim])**2) - self.lengths[i]**2
    return c

def linkage_C_fun(q, links, dim=3):
    """
    Compute Jacobian matrix of c at q.
    """    
    m = len(links)
    n = len(q)                                    
    C = np.zeros((m,n))
    for k, link in enumerate(links):
        for d in range(dim):
            C[k,link[0]*dim + d] += 2.0*(q[link[0]*dim + d] - q[link[1]*dim + d])
            C[k,link[1]*dim + d] += -2.0*(q[link[0]*dim + d] - q[link[1]*dim + d])
    return C
