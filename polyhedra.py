import numpy as np
import scipy as sp
import os

from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cm


#give coordinates of polyhedra centered at origin with unit edges and face centers, face inds

def check_face_index_order(polyf):
    """
    Make sure that the vertex indices are in clockwise order. 
    We assume the origin is inside the polyhedron.
    """
    
    # Get poly info.
    v, f_i, c = polyf()

    # Make sure poly is convex, or else we cannot guarentee clockwise order
    if is_convex(v, f_i) == False:
        print "WARNING: Polyhedron non-convex. Cannot guarantee correct ordering."
        return True

    bad_faces = []
    # Check for clockwise-ness
    for k, face in enumerate(f_i):
        if np.dot(v[face[1]],np.cross((v[face[0]] - v[face[1]]),(v[face[2]] - v[face[1]]))) < 0.0:
            print v[face[1]],np.cross((v[face[0]] - v[face[1]]),(v[face[2]] - v[face[1]]))
            bad_faces.append(k)
    if len(bad_faces) > 0:
            print "ERROR: these faces are not clockwise:", bad_faces
            return False
    
    return True

def is_convex(vs, f_i):
    """
    We assume that the polyhedron contains the origin and check for convexity."
    """

    # For each face, make sure all other vertices 
    # are on one side of the plane the face lays in.
    for face in f_i:
        #pos = 0
        #neg = 0
        pos = []
        neg = []

        for k,v in enumerate(vs):
            # Don't include vertices in current face
            if k in face:
                continue
            # Check if all other vertices lay on one side of the face
            if np.dot(v - vs[face[1]],np.cross((vs[face[0]] - vs[face[1]]),(vs[face[2]] - vs[face[1]]))) > 0.0:
                #pos += 1
                pos.append(k)
            else:
                #neg += 1
                neg.append(k)

        #if pos != 0 and neg != 0:
        if len(pos) != 0 and len(neg) != 0:
            print face, pos, neg
            return False

    return True



def get_face_centers(f,v):
    #N_f = f.shape[0]
    #N_s = f[0].shape[0]
    N_f = len(f)
    #N_s = len(f[0])
    centers = np.zeros((N_f,3))
    for k in range(N_f):
        fc = np.zeros((1,3))
        N_s = len(f[k])
        for j in range(N_s):
            #print j,k,N_f,N_s
            fc += v[f[k][j],:]
        fc /= N_s
        centers[k,:] = fc
    return centers
    

def plot_polyhedron(poly_fun):
    verts, face_inds, cents = poly_fun()

    ax = Axes3D(plt.figure())
    scale = np.abs(verts).max()*1.2
    ax.set_xlim(-scale,scale)
    ax.set_ylim(-scale,scale)
    ax.set_zlim(-scale,scale)
    for i in range(len(face_inds)):
        side = []
        for j in range(len(face_inds[i])):
            side.append([verts[face_inds[i][j],0],verts[face_inds[i][j],1],verts[face_inds[i][j],2]])

        tri = Poly3DCollection([side])
        color = colors.rgb2hex(sp.rand(3))
        tri.set_facecolor(color)
        tri.set_edgecolor('k')

        ax.add_collection3d(tri)

    plt.show()

def plot_vertex_labels(poly_fun,faces=False,centers=False):
    verts, face_inds, cents = poly_fun()

    ax = Axes3D(plt.figure())
    scale = np.abs(verts).max()*1.2
    ax.set_xlim(-scale,scale)
    ax.set_ylim(-scale,scale)
    ax.set_zlim(-scale,scale)
        
    ax.scatter(verts[:,0], verts[:,1], verts[:,2], c='r', marker='o',s=50)
    
    for i in range(len(verts)):
        ax.text3D(verts[i,0],verts[i,1],verts[i,2], " "+str(i),size="large")

    if faces == True:
        for i in range(len(face_inds)):
            side = []
            for j in range(len(face_inds[i])):
                side.append([verts[face_inds[i][j],0],verts[face_inds[i][j],1],verts[face_inds[i][j],2]])

            tri = Poly3DCollection([side])
            color = colors.rgb2hex(sp.rand(3))
            tri.set_facecolor(color)
            tri.set_edgecolor('k')

            ax.add_collection3d(tri)
    
    if centers == True:
        ax.scatter(cents[:,0], cents[:,1], cents[:,2], c='b', marker='o',s=50)
    
    
    plt.show()

def write_bg_input_file(poly_str):
    """
    Takes a string with a polyhedron's name and writes 
    the Building Game input file for that polyhedron
    """
    # Get poly embedding info
    #v, f, c = poly_fun()
    #v, f, c = getattr(polyhedra,poly_name)()
    v, f, c = globals()[poly_str]()
    
    # Compute polyhedron statistics
    V = len(v)
    F = len(f)
    E = V + F - 2

    # Find the types of each face
    species, f_types = get_face_type_info(f)
 
    # Compute the adjacency list
    adj_list = get_adj_list(v,f)

    # Compute dual list (faces adj to each vertex)
    dual = get_dual(v,f,adj_list)
    if dual == False:
        print "ERROR: could not compute dual. File was not written"
        return 

    # Write file
    filename = os.path.join(os.path.dirname(__file__),'data',poly_str + "_5_1.txt")
    
    try: 
        f_write = open(filename,'w')
    except:
        print "ERROR: Bad filename"
        return

    f_write.write(poly_str+'\n')
    f_write.write(str(F)+" "+str(E)+" "+str(V)+'\n')
    f_write.write(str(len(species)))
    for s in species:
        f_write.write(" "+str(s))
    f_write.write('\n')

    for j in range(len(adj_list)):
        f_write.write('1 '+str(f_types[j]))
        for a in adj_list[j]:
            f_write.write(' '+str(a))
        f_write.write('\n')

    for j in range(len(dual)):
        f_write.write(str(len(dual[j])))
        for b in dual[j]:
            f_write.write(' '+str(b))
        f_write.write('\n')

    f_write.close()
    
    return

     


def get_dual(v,f,adj_list):
    """
    For each vertex, make a list of faces that share this vertex
    """

    dual = []

    for vert in range(len(v)):
        v_dual = []
        
        for k,face in enumerate(f):
            #print vert, face
            if vert in face:
                v_dual.append(k)
        ordered_v_dual = order_v_dual(v_dual,adj_list)

        if ordered_v_dual != False:
            dual.append(ordered_v_dual)
        else:
            print "Dual error in vertex", vert
            print "Unordered dual:", v_dual
            print "Adjacency List:", adj_list
            return False
    return dual
        
def order_v_dual(y, adj_list):
    """
    Take the list of faces x adjacent to vertex vert 
    and order them in a clockwise fashion.
    """

    # Seed ordered list z with first element of the unordered list x
    x = y[:]

    if len(x) == 0:
        print "ERROR: empty dual"
        return False

    z = [x[0]]
    del x[0]

    # At each iteration determine the next clockwise face. 
    #Add it to z and remove it from x.
    while len(x) > 1:
        count = len(x)        
        #print x, z
        #print adj_list
        for k in range(-2, len(adj_list[z[-1]]) - 2):
            if adj_list[z[-1]][k] in y and adj_list[z[-1]][k + 1] in y:
                z.append(adj_list[z[-1]][k])
                x.remove(z[-1])
                break
        if count == len(x):
            print "ERROR: Removal failed in order_v_dual for dual",x
            return False

    # Append last remaining element
    #print "end",x,z
    z.append(x[0])

    return z
    
def get_face_type_info(f):
    """
    Take f and compute the different face species
    """
    # Get species
    species = list(set([len(face) for face in f]))

    # Get each faces species
    f_types = []
    for face in f:
        f_types.append(species.index(len(face)))

    return species, f_types
        


def get_adj_list(v,f):
    """
    For the vertex locations v and the list of each of the vertices 
    of each face f, create the face adjacency list. 
    """
    adj_list = []
    
    # For each face, find the adjacent ones in clockwise order
    for j in range(len(f)):
        f_adj = []
        for k in range(len(f[j])):
            # Get other face that shares the two vertices 
            f_adj.append(get_adj_face(f, j, f[j][k-1], f[j][k]))
        adj_list.append(f_adj)

    return adj_list

def get_adj_face(f, j, v1, v2):
    
    for k,face in enumerate(f):
        if k == j:
            continue
        if (v1 in face and v2 in face):
            return k

    print "Failed to find adjacent face."
    return False
    
    

def tetrahedron():
    verts = np.zeros((4,3))
    verts[0,:] = .5*np.array([1.0,	0.0,	-1/2**.5])
    verts[1,:] = .5*np.array([-1.0,	0.0,	-1/2**.5])
    verts[2,:] = .5*np.array([0.0,	1.0,	1/2**.5])
    verts[3,:] = .5*np.array([0.0,	-1.0,	1/2**.5])
    
    face_inds = [[0,	2,	1],
                 [0,	1,	3],
                 [0,	3,	2],
                 [1,	2,	3]]

    cents = get_face_centers(face_inds,verts)

    return verts, face_inds, cents

def cube():
    verts = np.zeros((8,3))
    verts[0,:] = .5*np.array([1.0,	1.0,	1.0])
    verts[1,:] = .5*np.array([1.0,	1.0,	-1.0])
    verts[2,:] = .5*np.array([1.0,	-1.0,	1.0])
    verts[3,:] = .5*np.array([1.0,	-1.0,	-1.0])
    verts[4,:] = .5*np.array([-1.0,	1.0,	1.0])
    verts[5,:] = .5*np.array([-1.0,	1.0,	-1.0])
    verts[6,:] = .5*np.array([-1.0,	-1.0,	1.0])
    verts[7,:] = .5*np.array([-1.0,	-1.0,	-1.0])

    face_inds = [[1,	5,	7,	3],
                 [0,	2,	6,	4],
                 [0,	1,	3,	2],
                 [2,	3,	7,	6],
                 [0,	4,	5,	1],
                 [4,	6,	7,	5]]
    
    cents = get_face_centers(face_inds,verts)

    return verts, face_inds, cents


def octahedron():
    verts = np.zeros((6,3))
    verts[0,:] = 2.0**-.5*np.array([1.0,	0.0,	0.0])
    verts[1,:] = 2.0**-.5*np.array([-1.0,	0.0,	0.0])
    verts[2,:] = 2.0**-.5*np.array([0.0,	1.0,	0.0])
    verts[3,:] = 2.0**-.5*np.array([0.0,	-1.0,	0.0])
    verts[4,:] = 2.0**-.5*np.array([0.0,	0.0,	1.0])
    verts[5,:] = 2.0**-.5*np.array([0.0,	0.0,	-1.0])
    
    face_inds = [[0,	3,	4],
                 [1,	4,	3],
                 [0,	4,	2],
                 [0,	5,	3],
                 [1,	3,	5],
                 [1,	2,	4],
                 [0,	2,	5],
                 [1,	5,	2]]
    
    cents = get_face_centers(face_inds,verts)

    return verts, face_inds, cents

def dodecahedron():
    phi = .5*(5**.5 + 1.0)
    verts = np.zeros((20,3))
    verts[0,:] = .5*phi*np.array([1.0,	1.0,	1.0])
    verts[1,:] = .5*phi*np.array([1.0,	1.0,	-1.0])
    verts[2,:] = .5*phi*np.array([1.0,	-1.0,	1.0])
    verts[3,:] = .5*phi*np.array([1.0,	-1.0,	-1.0])
    verts[4,:] = .5*phi*np.array([-1.0,	1.0,	1.0])
    verts[5,:] = .5*phi*np.array([-1.0,	1.0,	-1.0])
    verts[6,:] = .5*phi*np.array([-1.0,	-1.0,	1.0])
    verts[7,:] = .5*phi*np.array([-1.0,	-1.0,	-1.0])

    verts[8,:] = .5*phi*np.array([0.0,	1.0/phi,	phi])
    verts[9,:] = .5*phi*np.array([0.0,	1.0/phi,	-phi])
    verts[10,:] = .5*phi*np.array([0.0,	-1.0/phi,	phi])
    verts[11,:] = .5*phi*np.array([0.0,	-1.0/phi,	-phi])
    
    verts[12,:] = .5*phi*np.array([1.0/phi,	phi,	0.0])
    verts[13,:] = .5*phi*np.array([1.0/phi,	-phi,	0.0])
    verts[14,:] = .5*phi*np.array([-1.0/phi,	phi,	0.0])
    verts[15,:] = .5*phi*np.array([-1.0/phi,	-phi,	0.0])
    
    verts[16,:] = .5*phi*np.array([phi,	0.0,	1.0/phi])
    verts[17,:] = .5*phi*np.array([-phi,	0.0,	1.0/phi])
    verts[18,:] = .5*phi*np.array([phi,	0.0,	-1.0/phi])
    verts[19,:] = .5*phi*np.array([-phi,	0.0,	-1.0/phi])
    
    face_inds = [[12,	1, 	18,	16,	0], 
                 [16,	2, 	10,	8, 	0],
		 [16,	18,	3, 	13,	2],
		 [8,	4, 	14,	12,	0],
		 [8,	10,	6, 	17,	4],
		 [13,	15,	6, 	10,	2],
		 [11,	7, 	15,	13,	3],
		 [9,	11,	3, 	18,	1],
		 [12,	14,	5, 	9, 	1],
		 [17,	19,	5, 	14,	4],
		 [15,	7, 	19,	17,	6],
		 [19,	7, 	11,	9, 	5]]

    cents = get_face_centers(face_inds,verts)

    return verts, face_inds, cents


def icosahedron():
    phi = .5*(5**.5 + 1.0)
    verts = np.zeros((12,3))
    verts[0,:] = .5*np.array([0.0,	1.0,	phi])
    verts[1,:] = .5*np.array([0.0,	1.0,	-phi])
    verts[2,:] = .5*np.array([0.0,	-1.0,	phi])
    verts[3,:] = .5*np.array([0.0,	-1.0,	-phi])

    verts[4,:] = .5*np.array([1.0,	phi,	0.0])
    verts[5,:] = .5*np.array([1.0,	-phi,	0.0])
    verts[6,:] = .5*np.array([-1.0,	phi,	0.0])
    verts[7,:] = .5*np.array([-1.0,	-phi,	0.0])

    verts[8,:] = .5*np.array([phi,	0.0,	1.0])
    verts[9,:] = .5*np.array([phi,	0.0,	-1.0])
    verts[10,:] = .5*np.array([-phi,	0.0,	1.0])
    verts[11,:] = .5*np.array([-phi,	0.0,	-1.0])
    
    face_inds = [[0,	2,	10],
                 [0,	8,	2],
                 [2,	7,	10],
                 [0,	10,	6],
                 [0,	6,	4],
                 [0,	4,	8],
                 [4,	9,	8],
                 [5,	8,	9],
                 [2,	8,	5],
                 [2,	5,	7],
                 [3,	7,	5],
                 [3,	11,	7],
                 [7,	11,	10],
                 [6,	10,	11],
                 [1,	6,	11],
                 [1,	4,	6],
                 [1,	9,	4],
                 [1,	3,	9],
                 [3,	5,	9],
                 [1,	11,	3]]
    
    cents = get_face_centers(face_inds,verts)

    return verts, face_inds, cents

def truncated_tetrahedron():
    verts = np.zeros((12,3))
    verts[0,:] = 8.0**-0.5*np.array([3.0,	1.0,	1.0])
    verts[1,:] = 8.0**-0.5*np.array([1.0,	3.0,	1.0])
    verts[2,:] = 8.0**-0.5*np.array([1.0,	1.0,	3.0])

    verts[3,:] = 8.0**-0.5*np.array([-3.0,	-1.0,	1.0])
    verts[4,:] = 8.0**-0.5*np.array([-1.0,	-3.0,	1.0])
    verts[5,:] = 8.0**-0.5*np.array([-1.0,	-1.0,	3.0])

    verts[6,:] = 8.0**-0.5*np.array([-3.0,	1.0,	-1.0])
    verts[7,:] = 8.0**-0.5*np.array([-1.0,	3.0,	-1.0])
    verts[8,:] = 8.0**-0.5*np.array([-1.0,	1.0,	-3.0])

    verts[9,:] = 8.0**-0.5*np.array([3.0,	-1.0,	-1.0])
    verts[10,:] = 8.0**-0.5*np.array([1.0,	-3.0,	-1.0])
    verts[11,:] = 8.0**-0.5*np.array([1.0,	-1.0,	-3.0])

    
    face_inds = [[0,	2,	1],
                 [0,	1,	7,	8,	11,	9],
                 [0,	9,	10,	4,	5,	2],
                 [1,	2,	5,	3,	6,	7],
                 [6,	8,	7],
                 [9,	11,	10],
                 [3,	5,	4],
                 [3,	4,	10,	11,	8,	6]]
    
    cents = get_face_centers(face_inds,verts)

    return verts, face_inds, cents


def cuboctahedron():
    verts = np.zeros((12,3))

    verts[0,:] = 2.0**-0.5*np.array([1.0,	1.0,	0.0])
    verts[1,:] = 2.0**-0.5*np.array([1.0,	-1.0,	0.0])
    verts[2,:] = 2.0**-0.5*np.array([-1.0,	1.0,	0.0])
    verts[3,:] = 2.0**-0.5*np.array([-1.0,	-1.0,	0.0])

    verts[4,:] = 2.0**-0.5*np.array([1.0,	0.0,	1.0])
    verts[5,:] = 2.0**-0.5*np.array([1.0,	0.0,	-1.0])
    verts[6,:] = 2.0**-0.5*np.array([-1.0,	0.0,	1.0])
    verts[7,:] = 2.0**-0.5*np.array([-1.0,	0.0,	-1.0])

    verts[8,:] = 2.0**-0.5*np.array([0.0,	1.0,	1.0])
    verts[9,:] = 2.0**-0.5*np.array([0.0,	1.0,	-1.0])
    verts[10,:] = 2.0**-0.5*np.array([0.0,	-1.0,	1.0])
    verts[11,:] = 2.0**-0.5*np.array([0.0,	-1.0,	-1.0])

    
    face_inds = [[0,	4,	8],
                 [4,	10,	6,	8],
                 [0,	8,	2,	9],
                 [0,	5,	1,	4],
                 [1,	10,	4],
                 [3,	6,	10],
                 [2,	8,	6],
                 [2,	7,	9],
                 [0,	9,	5],
                 [1,	5,	11],
                 [1,	11,	3,	10],
                 [2,	6,	3,	7],
                 [5,	9,	7,	11],
                 [3,	11,	7]]
    
    cents = get_face_centers(face_inds,verts)

    return verts, face_inds, cents


def grid22():
    verts = np.zeros((12,3))

    verts[0,:] = np.array([0.0,		0.0,  	0.0])
    verts[1,:] = np.array([0.0,		1.0,	0.0])
    verts[2,:] = np.array([0.0,     	2.0,	0.0])
    verts[3,:] = np.array([1.0,		0.0,	0.0])
    verts[4,:] = np.array([1.0,		1.0,	0.0])
    verts[5,:] = np.array([1.0,		2.0,	0.0])
    verts[6,:] = np.array([2.0,		0.0,	0.0])
    verts[7,:] = np.array([2.0,		1.0,	0.0])
    verts[8,:] = np.array([2.0,		2.0,	0.0])
 
    
    face_inds = [[0,	1,	4,	3],
                 [1,	2,	5,	4],
                 [3,	4,	7,	6],
                 [4,	5,	8,	7]]
    
    cents = get_face_centers(face_inds,verts)

    return verts, face_inds, cents

def grid23():
    verts = np.zeros((12,3))

    verts[0,:] = np.array([0.0,		0.0,  	0.0])
    verts[1,:] = np.array([0.0,		1.0,	0.0])
    verts[2,:] = np.array([0.0,     	2.0,	0.0])
    verts[3,:] = np.array([1.0,		0.0,	0.0])
    verts[4,:] = np.array([1.0,		1.0,	0.0])
    verts[5,:] = np.array([1.0,		2.0,	0.0])
    verts[6,:] = np.array([2.0,		0.0,	0.0])
    verts[7,:] = np.array([2.0,		1.0,	0.0])
    verts[8,:] = np.array([2.0,		2.0,	0.0])
    verts[9,:] = np.array([3.0,		0.0,	0.0])
    verts[10,:] = np.array([3.0,	1.0,	0.0])
    verts[11,:] = np.array([3.0,	2.0,	0.0])

    
    face_inds = [[0,	1,	4,	3],
                 [1,	2,	5,	4],
                 [3,	4,	7,	6],
                 [4,	5,	8,	7],
                 [6,	7,	10,	9],
                 [7,	8,	11,	10]]
    
    cents = get_face_centers(face_inds,verts)

    return verts, face_inds, cents

def grid23b0():
    verts = np.zeros((12,3))

    verts[0,:] = np.array([0.0,		0.0,  			0.0])
    verts[1,:] = np.array([0.0,		1.0,			0.0])
    verts[2,:] = np.array([0.0,     	1.0+0.5*2.0**0.5,	0.5*2.0**0.5])
    verts[3,:] = np.array([1.0,		0.0,			0.0])
    verts[4,:] = np.array([1.0,		1.0,			0.0])
    verts[5,:] = np.array([1.0,		1.0+0.5*2.0**0.5,	0.5*2.0**0.5])
    verts[6,:] = np.array([2.0,		0.0,			0.0])
    verts[7,:] = np.array([2.0,		1.0,			0.0])
    verts[8,:] = np.array([2.0,		1.0+0.5*2.0**0.5,	0.5*2.0**0.5])
    verts[9,:] = np.array([3.0,		0.0,			0.0])
    verts[10,:] = np.array([3.0,	1.0,			0.0])
    verts[11,:] = np.array([3.0,	1.0+0.5*2.0**0.5,	0.5*2.0**0.5])

    
    face_inds = [[0,	1,	4,	3],
                 [1,	2,	5,	4],
                 [3,	4,	7,	6],
                 [4,	5,	8,	7],
                 [6,	7,	10,	9],
                 [7,	8,	11,	10]]
    
    cents = get_face_centers(face_inds,verts)

    return verts, face_inds, cents

def grid23b1():
    verts = np.zeros((12,3))

    verts[0,:] = np.array([0.0,		0.0,  	0.0])
    verts[1,:] = np.array([0.0,		1.0,	0.0])
    verts[2,:] = np.array([0.0,     	2.0,	0.0])
    verts[3,:] = np.array([1.0,		0.0,	0.0])
    verts[4,:] = np.array([1.0,		1.0,	0.0])
    verts[5,:] = np.array([1.0,		2.0,	0.0])
    verts[6,:] = np.array([1.0+0.5*2.0**0.5,		0.0,	0.5*2.0**0.5])
    verts[7,:] = np.array([1.0+0.5*2.0**0.5,		1.0,	0.5*2.0**0.5])
    verts[8,:] = np.array([1.0+0.5*2.0**0.5,		2.0,	0.5*2.0**0.5])
    verts[9,:] = np.array( [1.0+0.5*2.0**0.5,		0.0,	1.0+0.5*2.0**0.5])
    verts[10,:] = np.array([1.0+0.5*2.0**0.5,		1.0,	1.0+0.5*2.0**0.5])
    verts[11,:] = np.array([1.0+0.5*2.0**0.5,		2.0,	1.0+0.5*2.0**0.5]) 

    
    face_inds = [[0,	1,	4,	3],
                 [1,	2,	5,	4],
                 [3,	4,	7,	6],
                 [4,	5,	8,	7],
                 [6,	7,	10,	9],
                 [7,	8,	11,	10]]
    
    cents = get_face_centers(face_inds,verts)

    return verts, face_inds, cents


def grid13():
    verts = np.zeros((12,3))

    verts[0,:] = np.array([0.0,		0.0,  	0.0])
    verts[1,:] = np.array([0.0,		1.0,	0.0])
    verts[2,:] = np.array([1.0,		0.0,	0.0])
    verts[3,:] = np.array([1.0,		1.0,	0.0])
    verts[4,:] = np.array([2.0,		0.0,	0.0])
    verts[5,:] = np.array([2.0,		1.0,	0.0])
    verts[6,:] = np.array([3.0,		0.0,	0.0])
    verts[7,:] = np.array([3.0,		1.0,	0.0])
 
    
    face_inds = [[0,	1,	3,	2],
                 [2,	3,	5,	4],
                 [4,	5,	7,	6]]
    
    cents = get_face_centers(face_inds,verts)

    return verts, face_inds, cents



def truncated_cube():
    verts = np.zeros((24,3))
    xi = 2.0**.5 - 1.0

    verts[0,:] = 0.5/xi*np.array([xi,		1.0,	1.0])
    verts[1,:] = 0.5/xi*np.array([xi,		1.0,	-1.0])
    verts[2,:] = 0.5/xi*np.array([xi,		-1.0,	1.0])
    verts[3,:] = 0.5/xi*np.array([xi,		-1.0,	-1.0])
    verts[4,:] = 0.5/xi*np.array([-xi,		1.0,	1.0])
    verts[5,:] = 0.5/xi*np.array([-xi,		1.0,	-1.0])
    verts[6,:] = 0.5/xi*np.array([-xi,		-1.0,	1.0])
    verts[7,:] = 0.5/xi*np.array([-xi,		-1.0,	-1.0])

    verts[8,:] = 0.5/xi*np.array([1.0,		xi,	1.0])
    verts[9,:] = 0.5/xi*np.array([1.0,		xi,	-1.0])
    verts[10,:] = 0.5/xi*np.array([1.0,		-xi,	1.0])
    verts[11,:] = 0.5/xi*np.array([1.0,		-xi,	-1.0])
    verts[12,:] = 0.5/xi*np.array([-1.0,	xi,	1.0])
    verts[13,:] = 0.5/xi*np.array([-1.0,	xi,	-1.0])
    verts[14,:] = 0.5/xi*np.array([-1.0,	-xi,	1.0])
    verts[15,:] = 0.5/xi*np.array([-1.0,	-xi,	-1.0])

    verts[16,:] = 0.5/xi*np.array([1.0,		1.0,	xi])
    verts[17,:] = 0.5/xi*np.array([1.0,		1.0,	-xi])
    verts[18,:] = 0.5/xi*np.array([1.0,		-1.0,	xi])
    verts[19,:] = 0.5/xi*np.array([1.0,		-1.0,	-xi])
    verts[20,:] = 0.5/xi*np.array([-1.0,	1.0,	xi])
    verts[21,:] = 0.5/xi*np.array([-1.0,	1.0,	-xi])
    verts[22,:] = 0.5/xi*np.array([-1.0,	-1.0,	xi])
    verts[23,:] = 0.5/xi*np.array([-1.0,	-1.0,	-xi])
    
    
    face_inds = [[6,	22,	14],
                 [0,	8,	10,	2,	6,	14,	12,	4],
                 [2,	18,	19,	3,	7,	23,	22,	6],
                 [12,	14,	22,	23,	15,	13,	21,	20],
                 [4,	12,	20],
                 [0,	16,    	8],
                 [2,	10,	18],
                 [3,	19,	11],
                 [7,	15,	23],
                 [5,	21,	13],
                 [0,	4,	20,	21,	5,	1,	17,	16],
                 [8,	16,	17,	9,	11,	19,	18,	10],
                 [1,	5,	13,	15,	7,	3,	11,	9],
                 [1,	9,	17]]
    
    cents = get_face_centers(face_inds,verts)

    return verts, face_inds, cents


def rhombicuboctahedron():
    verts = np.zeros((24,3))

    verts[0,:] = 0.5*np.array([1.0,	1.0,	(1.0+2.0**0.5)])
    verts[1,:] = 0.5*np.array([1.0,	1.0,	-(1.0+2.0**0.5)])
    verts[2,:] = 0.5*np.array([1.0,	-1.0,	(1.0+2.0**0.5)])
    verts[3,:] = 0.5*np.array([1.0,	-1.0,	-(1.0+2.0**0.5)])
    verts[4,:] = 0.5*np.array([-1.0,	1.0,	(1.0+2.0**0.5)])
    verts[5,:] = 0.5*np.array([-1.0,	1.0,	-(1.0+2.0**0.5)])
    verts[6,:] = 0.5*np.array([-1.0,	-1.0,	(1.0+2.0**0.5)])
    verts[7,:] = 0.5*np.array([-1.0,	-1.0,	-(1.0+2.0**0.5)])

    verts[8,:] =  0.5*np.array([1.0,	(1.0+2.0**0.5),		1.0])
    verts[9,:] =  0.5*np.array([1.0,	(1.0+2.0**0.5),		-1.0])
    verts[10,:] = 0.5*np.array([1.0,	-(1.0+2.0**0.5),	1.0])
    verts[11,:] = 0.5*np.array([1.0,	-(1.0+2.0**0.5),	-1.0])
    verts[12,:] = 0.5*np.array([-1.0,	(1.0+2.0**0.5),		1.0])
    verts[13,:] = 0.5*np.array([-1.0,	(1.0+2.0**0.5),		-1.0])
    verts[14,:] = 0.5*np.array([-1.0,	-(1.0+2.0**0.5),	1.0])
    verts[15,:] = 0.5*np.array([-1.0,	-(1.0+2.0**0.5),	-1.0])

    verts[16,:] = 0.5*np.array([(1.0+2.0**0.5),		1.0,	1.0])
    verts[17,:] = 0.5*np.array([(1.0+2.0**0.5),		1.0,	-1.0])
    verts[18,:] = 0.5*np.array([(1.0+2.0**0.5),		-1.0,	1.0])
    verts[19,:] = 0.5*np.array([(1.0+2.0**0.5),		-1.0,	-1.0])
    verts[20,:] = 0.5*np.array([-(1.0+2.0**0.5),	1.0,	1.0])
    verts[21,:] = 0.5*np.array([-(1.0+2.0**0.5),	1.0,	-1.0])
    verts[22,:] = 0.5*np.array([-(1.0+2.0**0.5),	-1.0,	1.0])
    verts[23,:] = 0.5*np.array([-(1.0+2.0**0.5),	-1.0,	-1.0]) 
    
    
    face_inds = [[0, 2, 6, 4],
                 [0, 4, 12, 8],
                 [0, 8, 16],
                 [0, 16, 18, 2],
                 [2, 18, 10],
                 [2, 10, 14, 6],
                 [6, 14, 22],
                 [4, 6, 22, 20],
                 [4, 20, 12],
                 [8, 12, 13, 9],
                 [8, 9, 17, 16],
                 [16, 17, 19, 18],
                 [10, 18, 19, 11],
                 [10, 11, 15, 14],
                 [14, 15, 23, 22],
                 [20, 22, 23, 21],
                 [12, 20, 21, 13],
                 [1, 17, 9],
                 [3, 11, 19],
                 [7, 23, 15],
                 [5, 13, 21],
                 [1, 9, 13, 5],
                 [1, 3, 19, 17],
                 [3, 7, 15, 11],
                 [5, 21, 23, 7],
                 [1, 5, 7, 3]]

    cents = get_face_centers(face_inds,verts)

    return verts, face_inds, cents


def truncated_octahedron():
    verts = np.zeros((24,3))

    verts[0,:] = 2.0**-0.5*np.array([0.0,	1.0,	2.0])
    verts[1,:] = 2.0**-0.5*np.array([0.0,	1.0,	-2.0])
    verts[2,:] = 2.0**-0.5*np.array([0.0,	-1.0,	2.0])
    verts[3,:] = 2.0**-0.5*np.array([0.0,	-1.0,	-2.0])
    verts[4,:] = 2.0**-0.5*np.array([0.0,	2.0,	1.0])
    verts[5,:] = 2.0**-0.5*np.array([0.0,	2.0,	-1.0])
    verts[6,:] = 2.0**-0.5*np.array([0.0,	-2.0,	1.0])
    verts[7,:] = 2.0**-0.5*np.array([0.0,	-2.0,	-1.0])

    verts[8,:] =  2.0**-0.5*np.array([1.0,	0.0,	2.0])
    verts[9,:] =  2.0**-0.5*np.array([1.0,	0.0,	-2.0])
    verts[10,:] = 2.0**-0.5*np.array([-1.0,	0.0,	2.0])
    verts[11,:] = 2.0**-0.5*np.array([-1.0,	0.0,	-2.0])
    verts[12,:] = 2.0**-0.5*np.array([1.0,	2.0,	0.0])
    verts[13,:] = 2.0**-0.5*np.array([1.0,	-2.0,	0.0])
    verts[14,:] = 2.0**-0.5*np.array([-1.0,	2.0,	0.0])
    verts[15,:] = 2.0**-0.5*np.array([-1.0,	-2.0,	0.0])

    verts[16,:] = 2.0**-0.5*np.array([2.0,	0.0,	1.0])
    verts[17,:] = 2.0**-0.5*np.array([2.0,	0.0,	-1.0])
    verts[18,:] = 2.0**-0.5*np.array([-2.0,	0.0,	1.0])
    verts[19,:] = 2.0**-0.5*np.array([-2.0,	0.0,	-1.0])
    verts[20,:] = 2.0**-0.5*np.array([2.0,	1.0,	0.0])
    verts[21,:] = 2.0**-0.5*np.array([2.0,	-1.0,	0.0])
    verts[22,:] = 2.0**-0.5*np.array([-2.0,	1.0,	0.0])
    verts[23,:] = 2.0**-0.5*np.array([-2.0,	-1.0,	0.0])
    
    
    face_inds = [[0, 10, 18, 22, 14, 4],
                 [0, 8, 2, 10],
                 [18, 23, 19, 22],
                 [4, 14, 5, 12],
                 [2, 6, 15, 23, 18, 10],
                 [1, 5, 14, 22, 19, 11],
                 [0, 4, 12, 20, 16, 8],
                 [2, 8, 16, 21, 13, 6],
                 [3, 11, 19, 23, 15, 7],
                 [1, 9, 17, 20, 12, 5],
                 [6, 13, 7, 15],
                 [1, 11, 3, 9],
                 [16, 20, 17, 21],
                 [3, 7, 13, 21, 17, 9]]
    
    cents = get_face_centers(face_inds,verts)
    #cents = []
    return verts, face_inds, cents


def triakis_tetrahedron():
    verts = np.zeros((8,3))

    verts[0,:] = 8.0**-0.5*np.array([(5.0/3.0),	(5.0/3.0),	(5.0/3.0)])
                 
    verts[1,:] = 8.0**-0.5*np.array([1.0,	1.0,	-1.0])
    verts[2,:] = 8.0**-0.5*np.array([1.0,	-1.0,	1.0])
    verts[3,:] = 8.0**-0.5*np.array([-1.0,	1.0,	1.0])
                 
    verts[4,:] = 8.0**-0.5*np.array([-(5.0/3.0),	(5.0/3.0),	-(5.0/3.0)])
    verts[5,:] = 8.0**-0.5*np.array([(5.0/3.0),		-(5.0/3.0),	-(5.0/3.0)])
    verts[6,:] = 8.0**-0.5*np.array([-(5.0/3.0),	-(5.0/3.0),	(5.0/3.0)])
                 
    verts[7,:] = 8.0**-0.5*np.array([-1.0,	-1.0,	-1.0])
    
    
    face_inds = [[3, 6, 4],
                 [0, 3, 4],
                 [0, 6, 3],
                 [4, 6, 7],
                 [0, 4, 1],
                 [0, 2, 6],
                 [5, 7, 6],
                 [4, 7, 5],
                 [1, 4, 5],
                 [0, 1, 5],
                 [0, 5, 2],
                 [2, 5, 6]]
    
    cents = get_face_centers(face_inds,verts)
    #cents = []
    return verts, face_inds, cents


def tetrakis_hexahedron():
    verts = np.zeros((14,3))

    verts[0,:] = 3.0**-0.5*np.array([-1.0,	1.0,	1.0])              
    verts[1,:] = 4.0**-0.5*np.array([0.0,	0.0,	2.0])
    verts[2,:] = 4.0**-0.5*np.array([-2.0,	0.0,	0.0])
    verts[3,:] = 4.0**-0.5*np.array([0.0,	2.0,	0.0])    
    verts[4,:] = 3.0**-0.5*np.array([-1.0,	-1.0,	1.0])
    verts[5,:] = 3.0**-0.5*np.array([-1.0,	1.0,	-1.0])
    verts[6,:] = 3.0**-0.5*np.array([1.0,	1.0,	1.0])
    verts[7,:] = 3.0**-0.5*np.array([1.0,	-1.0,	1.0])
    verts[8,:] = 3.0**-0.5*np.array([-1.0,	-1.0,	-1.0])              
    verts[9,:] = 3.0**-0.5*np.array([1.0,	1.0,	-1.0])
    verts[10,:] = 4.0**-0.5*np.array([0.0,	-2.0,	0.0])
    verts[11,:] = 4.0**-0.5*np.array([0.0,	0.0,	-2.0])    
    verts[12,:] = 4.0**-0.5*np.array([2.0,	0.0,	0.0])
    verts[13,:] = 3.0**-0.5*np.array([1.0,	-1.0,	-1.0])
    
    
    face_inds = [[2, 4, 8],
                 [2, 8, 5],
                 [0, 2, 5],
                 [0, 4, 2],
                 [0, 1, 4],
                 [4, 10, 8],
                 [5, 8, 11],
                 [0, 5, 3],
                 [0, 6, 1],
                 [1, 7, 4],
                 [4, 7, 10],
                 [8, 10, 13],
                 [8, 13, 11],
                 [5, 11, 9],
                 [3, 5, 9],
                 [0, 3, 6],
                 [1, 6, 7],
                 [7, 13, 10],
                 [9, 11, 13],
                 [3, 9, 6],
                 [6, 12, 7],
                 [7, 12, 13],
                 [9, 13, 12],
                 [6, 9, 12]]
    
    cents = get_face_centers(face_inds,verts)
    #cents = []
    return verts, face_inds, cents



#def deltoidal_icositetrahedron():
#    verts = np.zeros((14,3))
#
#    verts[0,:] = 3.0**-0.5*np.array([-1.0,	1.0,	1.0])              
#    verts[1,:] = 4.0**-0.5*np.array([0.0,	0.0,	2.0])
#    verts[2,:] = 4.0**-0.5*np.array([-2.0,	0.0,	0.0])
#    verts[3,:] = 4.0**-0.5*np.array([0.0,	2.0,	0.0])    
#    verts[4,:] = 3.0**-0.5*np.array([-1.0,	-1.0,	1.0])
#    verts[5,:] = 3.0**-0.5*np.array([-1.0,	1.0,	-1.0])
#    verts[6,:] = 3.0**-0.5*np.array([1.0,	1.0,	1.0])
#    verts[7,:] = 3.0**-0.5*np.array([1.0,	-1.0,	1.0])
#    verts[8,:] = 3.0**-0.5*np.array([-1.0,	-1.0,	-1.0])              
#    verts[9,:] = 3.0**-0.5*np.array([1.0,	1.0,	-1.0])
#    verts[10,:] = 4.0**-0.5*np.array([0.0,	-2.0,	0.0])
#    verts[11,:] = 4.0**-0.5*np.array([0.0,	0.0,	-2.0])    
#    verts[12,:] = 4.0**-0.5*np.array([2.0,	0.0,	0.0])
#    verts[13,:] = 3.0**-0.5*np.array([1.0,	-1.0,	-1.0])
#    
#    
#    face_inds = [[2, 4, 8],
#                 [2, 8, 5],
#                 [0, 2, 5],
#                 [0, 4, 2],
#                 [0, 1, 4],
#                 [4, 10, 8],
#                 [5, 8, 11],
#                 [0, 5, 3],
#                 [0, 6, 1],
#                 [1, 7, 4],
#                 [4, 7, 10],
#                 [8, 10, 13],
#                 [8, 13, 11],
#                 [5, 11, 9],
#                 [3, 5, 9],
#                 [0, 3, 6],
#                 [1, 6, 7],
#                 [7, 13, 10],
#                 [9, 11, 13],
#                 [3, 9, 6],
#                 [6, 12, 7],
#                 [7, 12, 13],
#                 [9, 13, 12],
#                 [6, 9, 12]]
#    
#    cents = get_face_centers(face_inds,verts)
#    #cents = []
#    return verts, face_inds, cents



def icosidodecahedron():
    verts = np.zeros((30,3))

    verts[0,:] = 3.0**-0.5*np.array([-1.0,	1.0,	1.0])              
    verts[1,:] = 4.0**-0.5*np.array([0.0,	0.0,	2.0])
    verts[2,:] = 4.0**-0.5*np.array([-2.0,	0.0,	0.0])
    verts[3,:] = 4.0**-0.5*np.array([0.0,	2.0,	0.0])    
    verts[4,:] = 3.0**-0.5*np.array([-1.0,	-1.0,	1.0])
    verts[5,:] = 3.0**-0.5*np.array([-1.0,	1.0,	-1.0])
    verts[6,:] = 3.0**-0.5*np.array([1.0,	1.0,	1.0])
    verts[7,:] = 3.0**-0.5*np.array([1.0,	-1.0,	1.0])
    verts[8,:] = 3.0**-0.5*np.array([-1.0,	-1.0,	-1.0])              
    verts[9,:] = 3.0**-0.5*np.array([1.0,	1.0,	-1.0])
    verts[10,:] = 4.0**-0.5*np.array([0.0,	-2.0,	0.0])
    verts[11,:] = 4.0**-0.5*np.array([0.0,	0.0,	-2.0])    
    verts[12,:] = 4.0**-0.5*np.array([2.0,	0.0,	0.0])
    verts[13,:] = 3.0**-0.5*np.array([1.0,	-1.0,	-1.0])
    

    phi = (1.0 + 5.0**0.5)*0.5

    verts[0,:] =  np.array([0.0,0.0,+phi])
    verts[1,:] =  np.array([0.0,0.0,-phi])
    verts[2,:] =  np.array([0.0,+phi,0.0])
    verts[3,:] =  np.array([0.0,-phi,0.0])
    verts[4,:] =  np.array([+phi,0.0,0.0])
    verts[5,:] =  np.array([-phi,0.0,0.0])

    verts[6,:] =  np.array([+0.5,+0.5*phi,+0.5*(1.0+phi)])
    verts[7,:] =  np.array([+0.5,+0.5*phi,-0.5*(1.0+phi)])
    verts[8,:] =  np.array([+0.5,-0.5*phi,+0.5*(1.0+phi)])
    verts[9,:] =  np.array([+0.5,-0.5*phi,-0.5*(1.0+phi)])
    verts[10,:] = np.array([-0.5,+0.5*phi,+0.5*(1.0+phi)])
    verts[11,:] = np.array([-0.5,+0.5*phi,-0.5*(1.0+phi)])
    verts[12,:] = np.array([-0.5,-0.5*phi,+0.5*(1.0+phi)])
    verts[13,:] = np.array([-0.5,-0.5*phi,-0.5*(1.0+phi)])

    verts[14,:] = np.array([+0.5*phi,+0.5*(1.0+phi),+0.5])
    verts[15,:] = np.array([+0.5*phi,+0.5*(1.0+phi),-0.5])
    verts[16,:] = np.array([+0.5*phi,-0.5*(1.0+phi),+0.5])
    verts[17,:] = np.array([+0.5*phi,-0.5*(1.0+phi),-0.5])
    verts[18,:] = np.array([-0.5*phi,+0.5*(1.0+phi),+0.5])
    verts[19,:] = np.array([-0.5*phi,+0.5*(1.0+phi),-0.5])
    verts[20,:] = np.array([-0.5*phi,-0.5*(1.0+phi),+0.5])
    verts[21,:] = np.array([-0.5*phi,-0.5*(1.0+phi),-0.5])
                            
    verts[22,:] = np.array([+0.5*(1.0+phi),+0.5,+0.5*phi])
    verts[23,:] = np.array([+0.5*(1.0+phi),+0.5,-0.5*phi])
    verts[24,:] = np.array([+0.5*(1.0+phi),-0.5,+0.5*phi])
    verts[25,:] = np.array([+0.5*(1.0+phi),-0.5,-0.5*phi])
    verts[26,:] = np.array([-0.5*(1.0+phi),+0.5,+0.5*phi])
    verts[27,:] = np.array([-0.5*(1.0+phi),+0.5,-0.5*phi])
    verts[28,:] = np.array([-0.5*(1.0+phi),-0.5,+0.5*phi])
    verts[29,:] = np.array([-0.5*(1.0+phi),-0.5,-0.5*phi])


    



    face_inds = [[8,24,16],
                 [16,24,4,25,17],
                 [3,16,17],
                 [3,20,12,8,16],
                 [8,12,0],
                 [0,6,22,24,8],
                 [4,24,22],
                 [4,23,25],
                 [4,22,14,15,23],
                 [6,14,22],
                 [0,10,6],
                 [6,10,18,2,14],
                 [2,15,14],
                 [7,23,15],
                 [2,19,11,7,15],
                 [2,18,19],
                 [11,19,27],
                 [5,27,19,18,26],
                 [1,11,27,29,13],
                 [5,29,27],
                 [1,7,11],
                 [9,17,25],
                 [1,9,25,23,7],
                 [1,13,9],
                 [13,29,21],
                 [5,28,20,21,29],
                 [12,20,28],
                 [3,21,20],
                 [0,12,28,26,10],
                 [10,26,18],
                 [3,17,9,13,21],
                 [5,26,28]]
    
    cents = get_face_centers(face_inds,verts)
    
    return verts, face_inds, cents



def rhombic_dodecahedron():
    verts = np.zeros((14,3))

    verts[0,:]  = np.array([+1.0,	+1.0,	+1.0])              
    verts[1,:]  = np.array([+1.0,	+1.0,	-1.0])              
    verts[2,:]  = np.array([+1.0,	-1.0,	+1.0])              
    verts[3,:]  = np.array([+1.0,	-1.0,	-1.0])              
    verts[4,:]  = np.array([-1.0,	+1.0,	+1.0])              
    verts[5,:]  = np.array([-1.0,	+1.0,	-1.0])              
    verts[6,:]  = np.array([-1.0,	-1.0,	+1.0])              
    verts[7,:]  = np.array([-1.0,	-1.0,	-1.0])              
    verts[8,:]   = np.array([+2.0,	+0.0,	+0.0])              
    verts[9,:]   = np.array([-2.0,	+0.0,	+0.0])              
    verts[10,:]  = np.array([+0.0,	+2.0,	+0.0])              
    verts[11,:]  = np.array([+0.0,	-2.0,	+0.0])              
    verts[12,:]  = np.array([+0.0,	+0.0,	+2.0])              
    verts[13,:]  = np.array([+0.0,	+0.0,	-2.0])              
      
    face_inds = [[2, 8, 3, 11],
                 [1, 13, 3, 8],
                 [3, 13, 7, 11],
                 [5, 9, 7, 13],
                 [1, 10, 5, 13],
                 [5, 10, 4, 9],
                 [6, 11, 7, 9],
                 [4, 12, 6, 9],
                 [11, 6, 12, 2],
                 [0, 8, 2, 12],
                 [0, 10, 1, 8],
                 [0, 12, 4, 10]]
    
    cents = get_face_centers(face_inds,verts)
    #cents = []
    return verts, face_inds, cents




def triakis_octahedron():
    return poly_from_dual(truncated_cube)



def poly_from_dual(polyf):
    """
    Use poly data from polyf and return poly data of polyf's dual
    """
    #print 1
    dual_v, dual_f_i, dual_c = polyf()
    #print 2
    verts = get_face_centers(dual_f_i, dual_v)
    #print 'verts', verts
    dual_adj_list = get_adj_list(dual_v, dual_f_i)
    #print 'dual_adj_list', dual_adj_list 
    face_inds = get_dual(dual_v, dual_f_i, dual_adj_list)
    #print 'face_inds', face_inds
    cents = get_face_centers(face_inds,verts)
    #print 5
    return verts, face_inds, cents


def deltoidal_icositetrahedron():
    return poly_from_dual(rhombicuboctahedron)

def pentagonal_icositetrahedron():
    return poly_from_dual(snub_cube)

def rhombic_triacontrahedron():
    return poly_from_dual(icosidodecahedron)
    
def truncated_cuboctahedron():
    verts = np.zeros((48,3))

    A = 1.0
    B = 1.0 + 2.0**0.5
    C = 1.0 + 2.0*2.0**0.5

    verts[0,:] = np.array([+A,	+B,	+C])
    verts[1,:] = np.array([+A,	+B,	-C])
    verts[2,:] = np.array([+A,	-B,	+C])
    verts[3,:] = np.array([+A,	-B,	-C])
    verts[4,:] = np.array([-A,	+B,	+C])
    verts[5,:] = np.array([-A,	+B,	-C])
    verts[6,:] = np.array([-A,	-B,	+C])
    verts[7,:] = np.array([-A,	-B,	-C])
    verts[8,:] = np.array([+A,	+C,	+B])
    verts[9,:] = np.array([+A,	-C,	+B])
    verts[10,:] = np.array([+A,	+C,	-B])
    verts[11,:] = np.array([+A,	-C,	-B])
    verts[12,:] = np.array([-A,	+C,	+B])
    verts[13,:] = np.array([-A,	-C,	+B])
    verts[14,:] = np.array([-A,	+C,	-B])
    verts[15,:] = np.array([-A,	-C,	-B])
    verts[16,:] = np.array([+B,	+A,	+C])
    verts[17,:] = np.array([+B,	+A,	-C])
    verts[18,:] = np.array([-B,	+A,	+C])
    verts[19,:] = np.array([-B,	+A,	-C])
    verts[20,:] = np.array([+B,	-A,	+C])
    verts[21,:] = np.array([+B,	-A,	-C])
    verts[22,:] = np.array([-B,	-A,	+C])
    verts[23,:] = np.array([-B,	-A,	-C])
    verts[24,:] = np.array([+B,	+C,	+A])
    verts[25,:] = np.array([+B,	-C,	+A])
    verts[26,:] = np.array([-B,	+C,	+A])
    verts[27,:] = np.array([-B,	-C,	+A])
    verts[28,:] = np.array([+B,	+C,	-A])
    verts[29,:] = np.array([+B,	-C,	-A])
    verts[30,:] = np.array([-B,	+C,	-A])
    verts[31,:] = np.array([-B,	-C,	-A])
    verts[32,:] = np.array([+C,	+A,	+B])
    verts[33,:] = np.array([-C,	+A,	+B])
    verts[34,:] = np.array([+C,	+A,	-B])
    verts[35,:] = np.array([-C,	+A,	-B])
    verts[36,:] = np.array([+C,	-A,	+B])
    verts[37,:] = np.array([-C,	-A,	+B])
    verts[38,:] = np.array([+C,	-A,	-B])
    verts[39,:] = np.array([-C,	-A,	-B])
    verts[40,:] = np.array([+C,	+B,	+A])
    verts[41,:] = np.array([-C,	+B,	+A])
    verts[42,:] = np.array([+C,	-B,	+A])
    verts[43,:] = np.array([-C,	-B,	+A])
    verts[44,:] = np.array([+C,	+B,	-A])
    verts[45,:] = np.array([-C,	+B,	-A])
    verts[46,:] = np.array([+C,	-B,	-A])
    verts[47,:] = np.array([-C,	-B,	-A])
          
    face_inds = [[9, 25, 29, 11, 15, 31, 27, 13],
                 [8, 12, 26, 30, 14, 10, 28, 24],
                 [32, 40, 44, 34, 38, 46, 42, 36],
                 [33, 37, 43, 47, 39, 35, 45, 41],
                 [0, 16, 20, 2, 6, 22, 18, 4],
                 [1, 5, 19, 23, 7, 3, 21, 17],
                 [4, 18, 33, 41, 26, 12],
                 [46, 38, 21, 3, 11, 29],
                 [0, 8, 24, 40, 32, 16],
                 [7, 23, 39, 47, 31, 15],
                 [2, 20, 36, 42, 25, 9],
                 [5, 14, 30, 45, 35, 19],
                 [6, 13, 27, 43, 37, 22],
                 [1, 17, 34, 44, 28, 10],
                 [47, 43, 27, 31],
                 [40, 24, 28, 44],
                 [21, 38, 34, 17],
                 [18, 22, 37, 33],
                 [42, 46, 29, 25],
                 [26, 41, 45, 30],
                 [2, 9, 13, 6],
                 [1, 10, 14, 5],
                 [3, 7, 15, 11],
                 [0, 4, 12, 8],
                 [16, 32, 36, 20],
                 [19, 35, 39, 23]]

    
    cents = get_face_centers(face_inds,verts)
        
    return verts, face_inds, cents

    return

def truncated_dodecahedron():
    return

def truncated_icosahedron():
    return

def snub_cube():
    xi = ((17.0 + 3.0*33.0**0.5)**(1.0/3.0) - (-17.0 + 3.0*33.0**0.5)**(1.0/3.0) - 1.0)/3.0
    verts = np.zeros((24,3))


    c1 =(3.0*33.0**0.5 + 17.0)**(1.0/3.0)
    c2 =(3.0*33.0**0.5 - 17.0)**(1.0/3.0)
    c3 = (199.0 + 3.0*33.0**0.5)**(1.0/3.0)
    c4 = (199.0 - 3.0*33.0**0.5)**(1.0/3.0)

    C1 = ((4.0 - c1 + c2)/12.0)**0.5
    C2 = ((2.0 + c1 - c2)/12.0)**0.5
    C3 = ((4.0 + c3 + c4)/12.0)**0.5

    verts[0,:]  = np.array([+C1,	+C2,	-C3])              
    verts[1,:]  = np.array([+C1,	-C2,	+C3])              
    verts[2,:]  = np.array([-C1,	+C2,	+C3])              
    verts[3,:]  = np.array([-C1,	-C2,	-C3])              
    verts[4,:]  = np.array([+C2,	-C3,	+C1])              
    verts[5,:]  = np.array([-C2,	+C3,	+C1])              
    verts[6,:]  = np.array([+C2,	+C3,	-C1])              
    verts[7,:]  = np.array([-C2,	-C3,	-C1])              
    verts[8,:]  = np.array([-C3,	+C1,	+C2])              
    verts[9,:]  = np.array([+C3,	+C1,	-C2])              
    verts[10,:]  = np.array([+C3,	-C1,	+C2])              
    verts[11,:]  = np.array([-C3,	-C1,	-C2])              
    verts[12,:]  = np.array([+C2,	+C1,	+C3])              
    verts[13,:]  = np.array([-C2,	+C1,	-C3])              
    verts[14,:]  = np.array([+C2,	-C1,	-C3])              
    verts[15,:]  = np.array([-C2,	-C1,	+C3])              
    verts[16,:]  = np.array([+C3,	+C2,	+C1])              
    verts[17,:]  = np.array([-C3,	-C2,	+C1])              
    verts[18,:]  = np.array([-C3,	+C2,	-C1])              
    verts[19,:]  = np.array([+C3,	-C2,	-C1])              
    verts[20,:]  = np.array([+C1,	+C3,	+C2])              
    verts[21,:]  = np.array([+C1,	-C3,	-C2])              
    verts[22,:]  = np.array([-C1,	-C3,	+C2])              
    verts[23,:]  = np.array([-C1,	+C3,	-C2])              

    face_inds = [[1, 4, 22],
                 [4, 21, 7, 22],
                 [4, 19, 21],
                 [3, 7, 21],
                 [3, 21, 14],
                 [3, 14, 0, 13],
                 [0, 14, 9],
                 [9, 14, 19],
                 [9, 19, 10, 16],
                 [10, 12, 16],
                 [14, 21, 19],
                 [4, 10, 19],
                 [1, 10, 4],
                 [1, 12, 10],
                 [1, 15, 2, 12],
                 [1, 22, 15],
                 [15, 22, 17],
                 [7, 17, 22],
                 [7, 11, 17],
                 [2, 15, 8],
                 [15, 17, 8],
                 [8, 17, 11, 18],
                 [3, 11, 7],
                 [3, 13, 11],
                 [11, 13, 18],
                 [13, 23, 18],
                 [0, 23, 13],
                 [0, 6, 23],
                 [0, 9, 6],
                 [6, 9, 16],
                 [2, 8, 5],
                 [2, 5, 20],
                 [5, 8, 18],
                 [5, 18, 23],
                 [5, 23, 6, 20],
                 [6, 16, 20],
                 [2, 20, 12],
                 [20, 16, 12]]
    
    cents = get_face_centers(face_inds,verts)

    return verts, face_inds, cents



def disdyakis_dodecahedron():
    return poly_from_dual(truncated_cuboctahedron)

def triakis_icosahedron():
    return poly_from_dual(truncated_dodecahedron)

def pentakis_dodecahedron():
    return poly_from_dual(truncated_icosahedron)
