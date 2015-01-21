import numpy as np
import scipy as sp

from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cm

def get_plane_equation(V):
    """
    Given a 3x3 ndarray with each row a vertex in R^3, 
    compute the equation for the plane the triangle lies in.
    Plane: N \dot X + d = 0
    """
    N = np.cross(V[1,:] - V[0,:], V[2,:] - V[0,:])
    d = -np.dot(N, V[0,:])
    return N, d
         
#def planar_distance(v1, N, d):
#    """
#    Given a 3x3 ndarray V with each row a vertex in R^3,
#    and the equation for a plane in R^3. 
#    Determine if all three verticies are on the same side of the plane.
#    """
       
def check_same_side(V, N, d):
    """
    Given a 3x3 ndarray V with each row a vertex in R^3,
    and the equation for a plane in R^3. 
    Determine if all three verticies are on the same side of the plane.
    """
    ### AS IS ALLOWS FOR POINTS TO LIE IN PLANE
    c = np.dot(V, N) + d
    if max(c) <= 0.0:
        return True, c
    elif min(c) >= 0.0:
        return True, c
    else:
        return False, c

def get_vertex_projections(D, V):
    """
    Given planar intersection direction D and 3x3 triablge vertices V,
    return the projection of each projection onto the line. 
    """
    max_ind = np.argmax(abs(D))
    return V[:,max_ind]

def get_isolated_vertex(c):
    """
    Return a permutation of the indices index of c with the first entry 
    corresponding to the value with a different signed value from the other two and a  
    """
    if c[0]*c[1] >= 0.0:
        return [2, 1, 0]
    elif c[0]*c[2] >= 0.0:
        return [1, 2, 0]
    else:
        return [0, 1, 2]

def get_line_location(p0, p1, c0, c1):
    return p0 + (p1 - p0)*c0/(c0 - c1)

def check_interval_intersection(t1a, t1b, t2a, t2b):
    """
    Determine if the intervals (t1a, t1b) and (t2a, t12b) intersect.
    """
    if min(t1a, t1b) < t2a and max(t1a, t1b) > t2a:
        return True
    elif min(t1a, t1b) < t2b and max(t1a, t1b) > t2b:
        return True
    elif min(t2a, t2b) < t1b and max(t2a, t2b) > t1b:
        return True
    else:
        return False

def rescale_triangle(V, scale):
    """
    Take triangle V and rescale each point by a factor of scale 
    toward the triangle's center of mass.
    """
    C = np.tile(V.mean(axis=0), (3,1))
    return C + scale*(V - C)
    

def triangle_intersection(V1, V2, scaling=1.0):
    """
    Check for intersection between two triangles sitting in R^3.
    V1 and V2 are 3x3 ndarrays with each row a vertex.
    """

    # Rescale triangles to avoid edge or corner touching.
    V1 = rescale_triangle(V1, scaling)
    V2 = rescale_triangle(V2, scaling)

    ### From Moller 1997
    # Compute plane equation of triangle 2.
    N2, d2 = get_plane_equation(V2)

    # Reject if all points of triangle 1 are on same side.
    same_side_1, c1 = check_same_side(V1, N2, d2)
    if same_side_1 == True:
        return False

    # Compute plane equation of triangle 1.
    N1, d1 = get_plane_equation(V1)

    # Reject if all points of triangle 2 are on same side.
    same_side_2, c2 = check_same_side(V2, N1, d1)
    if same_side_2 == True:
        return False

    # Compute intersection line and project onto largest axis.
    D = np.cross(N1, N2)
    P1 = get_vertex_projections(D, V1)
    P2 = get_vertex_projections(D, V2)

    # Compute the intervals for each triangle.
    c_ind_1 = get_isolated_vertex(c1)
    c_ind_2 = get_isolated_vertex(c2)
    t1a = get_line_location(P1[c_ind_1[0]], 
                            P1[c_ind_1[1]],
                            c1[c_ind_1[0]], 
                            c1[c_ind_1[1]])
    t1b = get_line_location(P1[c_ind_1[0]], 
                            P1[c_ind_1[2]],
                            c1[c_ind_1[0]], 
                            c1[c_ind_1[2]])
    t2a = get_line_location(P2[c_ind_2[0]], 
                            P2[c_ind_2[1]],
                            c2[c_ind_2[0]], 
                            c2[c_ind_2[1]])
    t2b = get_line_location(P2[c_ind_2[0]], 
                            P2[c_ind_2[2]],
                            c2[c_ind_2[0]], 
                            c2[c_ind_2[2]])
                            
    # Intersect the intervals
    return check_interval_intersection(t1a, t1b, t2a, t2b)

def get_face(x, f):
    """
    Given a flattened configuaration vector x, and indices for a triangle's corners f,
    return a 3x3 ndarray with the triangles coordinates.
    """
    t = np.zeros((3,3))
    for k in range(3):
        t[k,:] = x[3*f[k]:3*f[k]+3]
    return t

    
def random_triangle():
    return np.random.uniform(size=(3,3))


def plot_triangles(V1, V2):
    ax = Axes3D(plt.figure())
    scale = 1.0
    ax.set_xlim(0.0,scale)
    ax.set_ylim(0.0,scale)
    ax.set_zlim(0.0,scale)

    T1 = Poly3DCollection([V1])
    color = colors.rgb2hex(sp.rand(3))
    T1.set_facecolor(color)
    T1.set_edgecolor('k')
    ax.add_collection(T1)
    T2 = Poly3DCollection([V2])
    color = colors.rgb2hex(sp.rand(3))
    T2.set_facecolor(color)
    T2.set_edgecolor('k')
    ax.add_collection(T2)
    
    plt.show()

def test_code():
    V1 = random_triangle()
    V2 = random_triangle()
    print triangle_intersection(V1, V2)
    plot_triangles(V1, V2)

#    for i in range(len(face_inds)):
#        side = []
#        for j in range(len(face_inds[i])):
#            side.append([verts[face_inds[i][j],0],verts[face_inds[i][j],1],verts[face_inds[i][j],2]])
#
#        tri = Poly3DCollection([side])
#        color = colors.rgb2hex(sp.rand(3))
#        tri.set_facecolor(color)
#        tri.set_edgecolor('k')
#
#        ax.add_collection3d(tri)
#
#    plt.show()


