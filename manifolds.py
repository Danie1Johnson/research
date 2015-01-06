"""
NAME_boundary(x) returns a numpy array with an entry for each boundary. 
Each entry negative for being within the boundary and positive for outside the boundary.

NAME_boundary_normal(x, b_num) Assuming x is on the b_num-th boundary (or very close to it),
return a unitvector (of the smae length as x) that is inwardly normal to the b_num-th boundary at x.
"""
import numpy as np

def get_manifold(manifold_name, kwargs={}):
    try:
        c_fun = globals()[manifold_name+"_c"]
        C_fun = globals()[manifold_name+"_C"]
        c = lambda x: c_fun(x,**kwargs)
        C = lambda x: C_fun(x,**kwargs)
        return c, C
    except KeyError, NameError:
        print "ERROR:", manifold_name, "not found."

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
            
