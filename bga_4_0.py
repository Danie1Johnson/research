import numpy as np
import numpy.linalg
import scipy as sp
import scipy.sparse
import scipy.sparse.linalg
import scipy.optimize
import copy
#import scipy.misc
import os.path
#import cPickle
#from math import floor,cos,acos,sin
import sys
import math
from time import time








#import numpy as np
import matplotlib.pyplot as plt
#from time import time
import matplotlib.animation as animation

import manifold_reflected_brownian_motion as mrbm

mrbm = reload(mrbm)





#from mpl_toolkits.mplot3d import Axes3D
#from mpl_toolkits.mplot3d.art3d import Poly3DCollection
#import matplotlib.pyplot as plt
#import pandas as pd
#from pandas import Series, DataFrame

import polyhedra as poly

bg_version = "5_1"

def get_filename(poly_name, prefix, suffix):
    """
    Return relative path filename. 
    """
    return os.path.join(os.path.dirname(__file__),'data',prefix + poly_name + suffix)


def get_poly(poly_name):
    """"
    Return data from BG input file for poly
    """
    f = open(get_filename(poly_name, '', '_'+bg_version+'.txt'), 'r')
    # Shape name 
    line = [x for x in f.readline().split()]
    if line[0] != poly_name:
        print "oh dear.",line[0],"is not actually the same as",poly_name
    # V E F
    [F,E,V] = [int(x) for x in f.readline().split()]
    # face types
    line_S = [int(x) for x in f.readline().split()]
    S = line_S[0]
    species = line_S[1:]
    # Adjacency list
    f_types = []
    adj_list = []    
    for k in range(F):
        line = [int(x) for x in f.readline().split()]
        f_types.append(line[1:1+line[0]])
        adj_list.append(line[1+line[0]:])
    dual = []
    for j in range(V):
        line = [int(x) for x in f.readline().split()]
        dual.append(line[1:])

    return V,E,F,S,species,f_types,adj_list,dual
    

def get_bg_ss(poly_name):
    """
    Read BG state space from file.
    """
    # Load file.
    f = open(get_filename(poly_name,'bg_' + bg_version + '_', '_arc.txt'), 'r')
    
    # Read vertices
    vert_line = f.readline().split()
    V = int(vert_line[-1])
    ints = []
    ids = []
    paths = []
    shell_int = []
    shell_paths = []
    for k in range(V):
        line = [int(x) for x in f.readline().split()]
        if k != line[0]:
            print "OH DEAR."
        ids.append(line[1])
        ints.append(line[5:])
        paths.append(line[3])
        shell_int.append(line[2])
        shell_paths.append(line[4])
    ints = np.array(ints) # Switch to np arrays earlier??
    paths = np.array(paths)
    ids = np.array(ids)
    shell_int = np.array(shell_int)
    shell_paths = np.array(shell_paths)

    # Read edges.
    edge_line = f.readline().split()
    E = int(edge_line[-1])
    edges = []
    shell_edge = []
    for j in range(E):
        line = [int(x) for x in f.readline().split()]
        edges.append(line[0:2])
        shell_edge.append(line[2])
    edges = np.array(edges)
    shell_edge = np.array(shell_edge)
    
    return ints, ids, paths, shell_int, shell_paths, edges, shell_edge

def get_bg_log(poly_name):
    """
    Read info about building game for poly.
    """
    # Load file.
    f = open(get_filename(poly_name,'bg_' + bg_version + '_', '_log.txt'), 'r')
    
    vertices_line = f.readline().split()
    v_counts = np.array([int(x) for x in f.readline().split()])
    if sum(v_counts) != int(vertices_line[-1]):
        print "ERROR: Differing intermediate count and sum", sum(v_counts), int(vertices_line[-1])

    shell_vertices_line = f.readline().split()
    s_v_counts = np.array([int(x) for x in f.readline().split()])
    if sum(s_v_counts) != int(shell_vertices_line[-1]):
        print "ERROR: Differing shell_intermediate count and sum", sum(s_v_counts), int(shell_vertices_line[-1])

    edges_line = f.readline().split()
    e_counts = np.array([int(x) for x in f.readline().split()])
    if sum(e_counts) != int(edges_line[-1]):
        print "ERROR: Differing edge count and sum", sum(e_counts), int(edges_line[-1])

    shell_edges_line = f.readline().split()
    s_e_counts = np.array([int(x) for x in f.readline().split()])
    if sum(s_e_counts) != int(shell_edges_line[-1]):
        print "ERROR: Differing shell edge count and sum", sum(s_e_counts), int(shell_edges_line[-1])

    time_line = f.readline().split()
    t_counts = np.array([int(x) for x in f.readline().split()])
    if sum(t_counts) != int(time_line[-1]):
        print "ERROR: Differing time count and sum", sum(t_counts), int(time_line[-1])

    pathways = int(f.readline().split()[-1])
    shell_pathways = int(f.readline().split()[-1])
        
    return v_counts, s_v_counts, e_counts, s_e_counts, t_counts, pathways, shell_pathways

def print_oeis_sequences(poly_name):
    """
    Take in poly_name ant print the corresponding 
    OEIS internal format input strings for both 
    general and shellable intermediates.
    """
    v, s_v = get_bg_log(poly_name)[:2]

    print_oeis_format(v, sc=False)
    print ''
    print_oeis_format(s_v, sc=True)
    

def print_oeis_format(v, sc=False):
    """
    print sequence v in OEIS format.
    """
    if sc == True:
        sc_str = ' simply-connected '
    else:
        sc_str = ' '

    print '%I'
    print '%S ',
    for k, x in enumerate(v):
        # Don't include the 0th term
        if k == 0:
            continue
        sys.stdout.write(str(x))
        if k != len(v) - 1:
            sys.stdout.write(',')
    print ''
    print '%N Number of'+sc_str+'one-sided n-ominoes on the',
    for j, w in enumerate(poly_name.split('_')):
        if j == len(poly_name.split('_')) - 1:
            print w+'.'
        else:
            print w
    print '%K nonn,fini,full'
    print '%O 1'
    print '%A _Daniel Johnson_'
    
    
def get_open_edges(ints,adj_list):
    """
    For each intermediate in ints, return the number of open edges it has.
    """
    Es = []
   
    for inter in ints:
        E = 0
        for k in range(len(inter)):
            if inter[k] == 1:
               for j in adj_list[k]:
                   if inter[j] == 0:
                       E += 1.0
        Es.append(E)

    return np.array(Es)

def get_closed_edges(ints,adj_list):
    """
    For each intermediate in ints, return the number of closed edges it has.
    """
    Es = []
    
    for inter in ints:
        E = 0
        for k in range(len(inter)):
            if inter[k] == 1:
               for j in adj_list[k]:
                   if inter[j] == 1:
                       E += 0.5 # Beacause of double counting

        Es.append(E)

    return np.array(Es)


def get_degeneracies(ints,edges,adj_list,Rs=None):
    """
    For each BGSS connection in edges, compute the forward and backward degeneracy numbers.
    """

    if Rs == None:
        Rs = generate_rotations(adj_list)

    Ss = []
    Ts = []

    for e in edges:
        chi_j, chi_k = get_chis(ints[e[0]],ints[e[1]],Rs)
        
        sf = 0
        sb = 0

        for i in range(len(chi_j)):
            if chi_j[i] == 0:
                chi_j_ei = copy.copy(chi_j)
                chi_j_ei[i] = 1
                if same_int(chi_j_ei,chi_k,Rs) == True:
                    sf += 1

            if chi_k[i] == 1:
                chi_k_ei = copy.copy(chi_k)
                chi_k_ei[i] = 0
                if same_int(chi_j,chi_k_ei,Rs) == True:
                    sb += 1

        Ss.append(sf)
        Ts.append(sb)

    return np.array(Ss), np.array(Ts)


def generate_rotations(adj_list):
    """
    Generate a list of index permutations corresponding to the 
    rotation group of adj_list's polyhedron.
    """
    rotations = []
    
    F = len(adj_list)
    
    for j in range(F):
        for k in range(len(adj_list[j])):
            ind = rotate_ind(j,k,adj_list)
            if ind[0] == -2:
                continue
            rotations.append(ind)

    return rotations

def rotate_ind(piv,r,adj_list):
    """
    Get index permutation coresponding to moving face 0 to piv 
    and rotating that face in adj_list r times.
    """
    f = 0
    F = len(adj_list)
    
    prev = []
    nxt = []
    ind = [-1 for x in range(F)]
    
    if piv >= F:
        print "bad piv!"
    
    # For reindexing, recenter at piv
    ind[0] = piv
    f += 1

    # Number of faces adjacent to piv (must be the same for [0] and [piv])
    N = len(adj_list[0]) 
    
    # If base faces different, return error
    if len(adj_list[0]) != len(adj_list[piv]):
        return [-2]
    
    # Fill in indices for faces adjacent to piv
    for j in range(N):
        ind[adj_list[0][j]] = adj_list[piv][(r+j-1+N)%N]
        prev.append(adj_list[0][j])
        f += 1

    while f < F:
        if len(prev) == 0:
            print "ERROR, inf-loop"
            break
        
        for k in range(len(prev)):
            Nk = len(adj_list[prev[k]])

            for h in range(Nk):
                # Find already indexed face
                i = adj_list[prev[k]][h]
                if ind[i] != -1:
                    break
                
            for q in range(Nk):
                if adj_list[ind[prev[k]]][q] == ind[i]:
                    break

            for g in range(Nk):
                map_from = adj_list[prev[k]][(h+g)%Nk]
                map_to = adj_list[ind[prev[k]]][(q+g)%Nk]

                # Check if already indexed 
                if ind[map_from] != -1:
                    if ind[map_from] != map_to:
                        return [-2]

                else:
                    # If not, add index
                    ind[map_from] = map_to
                    nxt.append(map_from)
                    f += 1
              
        prev = []
        prev = copy.copy(nxt)
        nxt = []
    return ind

def get_chis(x1,x2,Rs):
    """
    Purmute the entries of x2 (according to a polyhedral rotation) such that 
    x1 is a sub int of x2.
    """
    if sum(x1) + 1 != sum(x2):
        print 'ERROR: non-coresponding x1 and x2, wrong number of faces.'

    chi_j = x1

    if sub_int(chi_j, x2) == True:
        chi_k = x2
        return chi_j, chi_k
    else:
        for R in Rs:
            chi_k = np.array([x2[R[k]] for k in range(len(R))])
            if sub_int(chi_j,chi_k) == True:
                return chi_j, chi_k

    print 'ERROR: non-coresponding x1 and x2 after search.'
    return

def sub_int(x1,x2):
    """
    Check if x1 == x2 in all entries except for exactly one.
    """
    if numpy.linalg.norm(x1 + x2) == (4*sum(x1)+1)**0.5:
        return True
    return False


def same_int(x1,x2,Rs):
    """
    Check if x1 is a polyhedral rotation of x2.
    """
    for R in Rs:
        if numpy.linalg.norm(np.array([x1[R[k]] for k in range(len(R))]) - x2) == 0:
            return True

    return False

def get_rs(ints,adj_list,Rs=None):
    """
    Compute the order of the rotation group for each intermediate.
    """
    if Rs == None:
        Rs = generate_rotations(adj_list)
    rs = []

    for inter in ints:
        r = 0
        for R in Rs:
            if numpy.linalg.norm(inter - np.array([inter[R[k]] for k in range(len(inter))])) == 0:
                r += 1
        rs.append(r)

    return np.array(rs)




def get_Q(BG_cons, betas=1.0, alpha=1.0):
    """
    
    """


    Q_dat = np.hstack((BG_cons['Q_jk'],
                       BG_cons['Q_kj'], 
                       -BG_ints['z_j']))

    Q_j = np.hstack((BG_cons['x_j'],
                     BG_cons['x_k'],
                     np.arange(N)))

    Q_k = np.hstack((BG_cons['x_k'],
                     BG_cons['x_j'],
                     np.arange(N)))
    
    #Q_dat = np.hstack((BG_cons[BG_cons['x_j'] != N - 1]['Q_jk'],
    #                   BG_cons[BG_cons['x_k'] != N - 1]['Q_kj'], 
    #                   -BG_ints['z_j'][:N-1].values))
    #
    #Q_j = np.hstack((BG_cons[BG_cons['x_j'] != N - 1]['x_j'],
    #                 BG_cons[BG_cons['x_k'] != N - 1]['x_k'],
    #                 np.arange(N - 1)))
    #
    #Q_k = np.hstack((BG_cons[BG_cons['x_j'] != N - 1]['x_k'],
    #                 BG_cons[BG_cons['x_k'] != N - 1]['x_j'],
    #                 np.arange(N - 1)))

    return scipy.sparse.coo_matrix((Q_dat,(Q_j,Q_k)), shape=(N,N)).tocsc()


def get_taus(Q, A=None):
    """
    """

    One_N = scipy.sparse.coo_matrix((np.array([1]),(np.array([N-1]),np.array([0]))), shape=(N,1))
    Ones_mN = scipy.sparse.coo_matrix((np.ones(N-1),(np.arange(N-1),np.zeros(N-1))), shape=(N,1))

    # Identity matrix excpt with M_j,j = 0 for if j not in A
    Ident_A = scipy.sparse.coo_matrix((np.array([1]),(np.array([N-1]),np.array([N-1]))), shape=(N,N))

    # Identity matrix excpt with M_j,j = 0 for if j in A
    Ident_Ac = scipy.sparse.coo_matrix((np.ones(N-1),(np.arange(N-1),np.arange(N-1))), shape=(N,N))

    # Column matrix with the jth row equal to one if j not in A
    Ones_Ac = scipy.sparse.coo_matrix((np.ones(N-1),(np.arange(N-1),np.zeros(N-1))), shape=(N,1))


    #Ident_A = Ident_A.tocsr()
    #Ident_Ac = Ident_Ac.tocsr()
    #Ones_Ac = Ones_Ac.tocsr()
    #Q = Q.tocsr()

    Ident_A = Ident_A.tocsc()
    Ident_Ac = Ident_Ac.tocsc()
    Ones_Ac = Ones_Ac.tocsc()
    Q = Q.tocsc()

    return scipy.sparse.linalg.spsolve(Ident_A - Ident_Ac*Q, Ones_Ac) 

def get_us(Q, t_start, t_stop, t_num, A=None):
    """
    """
    N = Q.shape[0]

    u0 = np.zeros(N)
    u0[N-1] = 1.0

    #One_N = scipy.sparse.coo_matrix((np.array([1]),(np.array([N-1]),np.array([0]))), shape=(N,1))
    #Ones_mN = scipy.sparse.coo_matrix((np.ones(N-1),(np.arange(N-1),np.zeros(N-1))), shape=(N,1))

    # Identity matrix excpt with M_j,j = 0 for if j not in A
    #Ident_A = scipy.sparse.coo_matrix((np.array([1]),(np.array([N-1]),np.array([N-1]))), shape=(N,N))

    # Identity matrix excpt with M_j,j = 0 for if j in A
    Ident_Ac = scipy.sparse.coo_matrix((np.ones(N-1),(np.arange(N-1),np.arange(N-1))), shape=(N,N))

    # Column matrix with the jth row equal to one if j not in A
    #Ones_Ac = scipy.sparse.coo_matrix((np.ones(N-1),(np.arange(N-1),np.zeros(N-1))), shape=(N,1))


    #Ident_A = Ident_A.tocsr()
    #Ident_Ac = Ident_Ac.tocsr()
    #Ones_Ac = Ones_Ac.tocsr()
    #Q = Q.tocsr()

    #Ident_A = Ident_A.tocsc()
    Ident_Ac = Ident_Ac.tocsc()
    #Ones_Ac = Ones_Ac.tocsc()
    Q = Q.tocsc()


    #u = scipy.sparse.linalg.expm_multiply(Ident_Ac*(Q - scipy.sparse.dia_matrix((Q.diagonal()[scipy.newaxis, :], [0]), shape=(N, N))), 
    #                                      u0, 
    #                                      start=t_start, 
    #                                      stop=t_stop, 
    #                                      num=t_num, 
    #                                      endpoint=True)
    
    u = scipy.sparse.linalg.expm_multiply(Ident_Ac*Q, 
                                          u0, 
                                          start=t_start, 
                                          stop=t_stop, 
                                          num=t_num, 
                                          endpoint=True)

    du = Ident_Ac*Q*u.T 


    #return scipy.sparse.linalg.spsolve(Ident_A - Ident_Ac*Q, Ones_Ac) 
    return u, du

def fit_gamma_dist(t,u):
    """
    Take a np arrays of positive time points t and function values u 
    and find parameters for the Gamm distribution that best fit (t,u)
    """

    beta = get_tail_rate(t, u)
    
    ## Function to minimize.
    #f = lambda x: numpy.linalg.norm(u - gamma_pdf(t, x[0], x[1]))
    #
    ## Initial Guess
    #x0 = np.array([1.0, 1.0])
    
    # Function to minimize.
    f = lambda x: numpy.linalg.norm(u[1:] - tau_pdf_est(t[1:], x[0], x[1], x[2]))
    
    # Initial Guess
    x0 = np.array([1.0, 0.0015, 1.0])

    res = scipy.optimize.minimize(f, x0)

    return res

def gamma_pdf(t, alpha, beta):
    """
    Return the gamma distribution PDF for each point in t.
    """
    
    return beta**alpha*t**(alpha - 1.0)*np.exp(-beta*t)/math.gamma(alpha)

def tau_pdf_est(t, gamma, omega, zeta):
    """
    Return the estimated PDF for each point in t.
    """

    return zeta*np.exp(-gamma/t - omega*t)


def get_tail_rate(t, u, interval_num=10):
    """
    Estimate the exponential decay parameter for the tail of u.
    """

    return -((np.log(u[-1]) - np.log(u[-interval_num]))/(t[-1] - t[-interval_num]))
    


def load_bg_int(poly_name, int_num):
    """
    Load information about specified polyhedral intermediate. 
    Return information for corresponding linkage.
    """

    try: 
        poly_info = getattr(poly, poly_name)
    except AttributeError:
        print "ERROR:", poly_name, "not found in polyhedra.py"
        raise

    verts, face_inds, cents = poly_info()
    V, E, F, S, species, f_types, adj_list, dual = get_poly(poly_name)
    ints, ids, paths, shell_int, shell_paths, edges, shell_edge = get_bg_ss(poly_name)

    try:
        int_faces = ints[int_num]
    except IndexError:
        print "ERROR:", poly_name, "does not have an intermediate", int_num
        raise

    # Reindex vertices in the intermediate.
    v_in_int = np.array([False for k in range(V)])
    for k in range(F):
        if int_faces[k] != 0:
            v_in_int[face_inds[k][0]] = True
            v_in_int[face_inds[k][1]] = True
            v_in_int[face_inds[k][2]] = True

    # Create table for translating between new and old vertex indexing.
    old_v_ind_from_new = []
    new_v_ind_from_old = []

    for v in range(V):
        if v_in_int[v] == True:
            new_v_ind_from_old.append(len(old_v_ind_from_new)) 
            old_v_ind_from_new.append(v)
        else:
            new_v_ind_from_old.append(-1)

    # Create faces, masses, links, and lengths 
    N = len(old_v_ind_from_new)
    dim = 3
    face_inds_new = [[new_v_ind_from_old[face_inds[j][k]] for k in range(3)] for j in range(V) if v_in_int[j]]
    q0 = np.array([])
    for k in range(N):
        for i in range(dim):
            q0 = np.hstack((q0,verts[old_v_ind_from_new[k],i]))
    masses = np.ones(N)
    
    links = []
    faces = []
    for f in range(F):
        if int_faces[f] != 0:
            new_face = []
            for j in range(len(face_inds[f])):
                a = face_inds[f][j-1]
                b = face_inds[f][j]
                a_new = new_v_ind_from_old[a]
                b_new = new_v_ind_from_old[b]
                mx = max(a_new, b_new)
                mn = min(a_new, b_new)
                #print hi, links, mn, mx
                if ([mn, mx] in links) == False:
                    links.append([mn, mx])
                new_face.append(b_new)
            faces.append(new_face)

    lengths = np.array([numpy.linalg.norm(verts[old_v_ind_from_new[links[k][0]],:] - verts[old_v_ind_from_new[links[k][1]],:]) for k in range(len(links))])


    return N, dim, q0, masses, links, lengths, faces


#def load_bg_int(poly_name, int_num):
#    """
#    Load information about specified polyhedral intermediate. 
#    Return information for corresponding linkage.
#    """
#
#    try: 
#        poly_info = getattr(poly, poly_name)
#    except AttributeError:
#        print "ERROR:", poly_name, "not found in polyhedra.py"
#        raise
#
#    verts, face_inds, cents = poly_info()
#    V, E, F, S, species, f_types, adj_list, dual = get_poly(poly_name)
#    ints, ids, paths, shell_int, shell_paths, edges, shell_edge = get_bg_ss(poly_name)
#
#    try:
#        int_faces = ints[int_num]
#    except IndexError:
#        print "ERROR:", poly_name, "does not have an intermediate", int_num
#        raise
#
#    # Reindex vertices in the intermediate.
#    v_in_int = np.array([False for k in range(V)])
#    for k in range(F):
#        if int_faces[k] != 0:
#            v_in_int[face_inds[k][0]] = True
#            v_in_int[face_inds[k][1]] = True
#            v_in_int[face_inds[k][2]] = True
#
#    # Create table for translating between new and old vertex indexing.
#    old_v_ind_from_new = []
#    new_v_ind_from_old = []
#
#    for v in range(V):
#        if v_in_int[v] == True:
#            new_v_ind_from_old.append(len(old_v_ind_from_new)) 
#            old_v_ind_from_new.append(v)
#        else:
#            new_v_ind_from_old.append(-1)
#
#    # Create faces, masses, links, and lengths 
#    N = len(old_v_ind_from_new)
#    dim = 3
#    face_inds_new = [[new_v_ind_from_old[face_inds[j][k]] for k in range(3)] for j in range(V) if v_in_int[j]]
#    q0 = np.array([])
#    for k in range(N):
#        for i in range(dim):
#            q0 = np.hstack((q0,verts[old_v_ind_from_new[k],i]))
#    masses = np.ones(N)
#    
#    links = []
#    faces = []
#    for f in range(F):
#        if int_faces[f] != 0:
#            new_face = []
#            for j in range(len(face_inds[f])):
#                a = face_inds[f][j-1]
#                b = face_inds[f][j]
#                a_new = new_v_ind_from_old[a]
#                b_new = new_v_ind_from_old[b]
#                mx = max(a_new, b_new)
#                mn = min(a_new, b_new)
#                #print hi, links, mn, mx
#                if ([mn, mx] in links) == False:
#                    links.append([mn, mx])
#                new_face.append(b_new)
#            faces.append(new_face)
#
#    lengths = np.array([numpy.linalg.norm(verts[old_v_ind_from_new[links[k][0]],:] - verts[old_v_ind_from_new[links[k][1]],:]) for k in range(len(links))])
#
#
#    return N, dim, q0, masses, links, lengths, faces
#

def face_position(bg_int, face_num, faces, dim=3):
    """
    Return the current and last positions in the desired dimensions for viewing.
    """
    x_list = []
    y_list = []     
     
    ## Right now, only views in x and y dimensions
    for v in faces[face_num]:
        x_list.append(bg_int.x[dim*v])
        y_list.append(bg_int.x[dim*v + 1])
    x_list.append(bg_int.x[dim*faces[face_num][0]])
    y_list.append(bg_int.x[dim*faces[face_num][0] + 1])
    x_list = np.array(x_list)
    y_list = np.array(y_list)
    return (x_list, y_list)


def init_aux(face_lines):
    """
    initialize animation
    """
    for k, f in enumerate(face_lines):
        f.set_data([], [])
    return face_lines

def animate_aux(i, bg_int, faces, face_lines):
    """
    perform animation step
    """
    bg_int.sample()
    for k, f in enumerate(face_lines):
        f.set_data(face_position(bg_int, k, faces))
    return face_lines

def bg_animation(bg_int, faces, save_animation=False, L=1.0):
    """
    """
    fig = plt.figure()
    ax = fig.add_subplot(111, 
                         aspect='equal', 
                         autoscale_on=False,
                         xlim=(-2.0, 2.0), 
                         ylim=(-2.0, 2.0))
    ax.grid()

    face_lines = ()
    for f in faces:
        temp_line, = ax.plot([], [], 'o-', lw=2)
        face_lines += (temp_line,)

    animate = lambda x: animate_aux(x, bg_int, faces, face_lines)
    init = lambda: init_aux(face_lines)

    # choose the interval based on dt and the time to animate one step
    t0 = time()
    animate(0)
    t1 = time()
    interval = 1000 * (L * bg_int.h) - (t1 - t0)

    ani = animation.FuncAnimation(fig, 
                                  animate, 
                                  frames=300,
                                  interval=interval, 
                                  blit=True, 
                                  init_func=init)

    if save_animation == True:
        # save the animation as an mp4.  This requires ffmpeg or mencoder to be
        # installed.  The extra_args ensure that the x264 codec is used, so that
        # the video can be embedded in html5.  You may need to adjust this for
        # your system: for more information, see
        # http://matplotlib.sourceforge.net/api/animation_api.html
        ani.save('triangular_linkage_diffusion.mp4', fps=3, extra_args=['-vcodec', 'libx264'])

    plt.show()

def find_dihedrals(f1, f2, x, faces):
    """
    Take the index of two faces of a BG intermediate and compute the dihedral angle between them
    """
    n1 = triangle_normal(f1, x, faces)
    n2 = triangle_normal(f2, x, faces)
    return np.arccos(np.dot(n1, n2))

def triangle_normal(f, x, faces):
    """
    Find normal vector to triangle. 
    """
    v1 = x[3*faces[f][0]:3*faces[f][0]+3] - x[3*faces[f][1]:3*faces[f][1]+3] 
    v2 = x[3*faces[f][2]:3*faces[f][2]+3] - x[3*faces[f][1]:3*faces[f][1]+3] 
    n = np.cross(v1, v2)
    return n/numpy.linalg.norm(n)

#
##poly_name = 'tetrahedron'
##poly_name = 'cube'
#poly_name = 'octahedron'
#
##poly_name = 'truncated_tetrahedron' ###########ERROR: non-coresponding x1 and x2, wrong number of faces.#####ERROR: non-coresponding x1 and x2 after search.
#
##poly_name = 'dodecahedron'
##poly_name = 'triakis_tetrahedron'
##poly_name = 'rhombic_dodecahedron'
##poly_name = 'cuboctahedron' ############
##poly_name = 'truncated_cube' ######
##poly_name = 'truncated_octahedron' #####
##poly_name = 'icosahedron'
##poly_name = 'triakis_octahedron'
##poly_name = 'tetrakis_hexahedron'
##poly_name = 'deltoidal_icositetrahedron'
##poly_name = 'pentagonal_icositetrahedron'
##poly_name = 'rhombicuboctahedron'
##poly_name = 'truncated_cuboctahedron'
##poly_name = 'rhombic_triacontahedron'
##poly_name = 'icosidodecahedron' ## LOG ONLY
#
#V, E, F, S, species, f_types, adj_list, dual = get_poly(poly_name)
#v_counts, s_v_counts, e_counts, s_e_counts, t_counts, pathways, shell_pathways = get_bg_log(poly_name)
#ints, ids, paths, shell_int, shell_paths, edges, shell_edge = get_bg_ss(poly_name)
#
#Rs = generate_rotations(adj_list)
#
#
#alpha = 3.0
#beta = 1.0
#
##print_oeis_sequences(poly_name)
#
#BG_ints = DataFrame(ids)
#BG_ints.columns = ['ids']
#BG_ints['shellable'] = shell_int
#BG_ints['paths'] = paths
#BG_ints['shell paths'] = shell_paths
#BG_ints['faces'] = ints.sum(axis=1)
#BG_ints['open edges'] = get_open_edges(ints, adj_list)
#BG_ints['closed edges'] = get_closed_edges(ints, adj_list)
#BG_ints['E_j'] = BG_ints['open edges']
#BG_ints['r_j'] = get_rs(ints, adj_list)
#BG_ints['pi_j'] = np.exp(-beta*BG_ints['E_j'])/BG_ints['r_j']
#
#z = sum(BG_ints['pi_j'])
#
#BG_ints['pi_j'] /= z
#
##for k in range(len(ints[0])):
##    BG_ints[str(k)] = ints[:,k]
#
#
#
#
#
##BG_ints.columns = [['Stats' for x in range(4)] + ['Intermediate' for y in range(len(ints[0]))],
##                   ['ids', 'shellable', 'paths', 'shell paths'] + [str(z) for z in range(len(ints[0]))]]
#
#    
##BG_ints.columns = [['Stats' for x in range(4)] + ['Intermediate' for y in range(len(ints[0]))],
##                   ['ids', 'shellable', 'paths', 'shell paths'] + [str(z) for z in range(len(ints[0]))]]
#
#
##betas = np.arange(0.0,2.0,0.02)
##alphas = np.arange(0.0,5.0,0.1)
##taus = []
##nus = []
#N = len(BG_ints)
#
#
#BG_cons = DataFrame(edges[:,0])
#BG_cons.columns = ['x_j']
#BG_cons['x_k'] = edges[:,1]
#BG_cons['shell edge'] = shell_edge
#BG_cons['S_jk'], BG_cons['S_kj'] =  get_degeneracies(ints, edges, adj_list)
#
##for beta in betas:    
##for alpha in alphas:    
#
#BG_cons['E_jk'] = np.vstack((BG_ints['E_j'][BG_cons['x_j']].values,
#                             BG_ints['E_j'][BG_cons['x_k']].values)).max(axis=0) + alpha
#
#
#BG_cons['B_jk'] = BG_cons['E_jk'].values - BG_ints['E_j'][BG_cons['x_j']].values
#BG_cons['B_kj'] = BG_cons['E_jk'].values - BG_ints['E_j'][BG_cons['x_k']].values
#
#
#BG_cons['Q_jk'] = BG_cons['S_jk'].values*np.exp(-beta*BG_cons['B_jk'].values)
#BG_cons['Q_kj'] = BG_cons['S_kj'].values*np.exp(-beta*BG_cons['B_kj'].values)
#
#BG_ints['z_j'] = np.array([BG_cons[BG_cons['x_j'] == j]['Q_jk'].sum() 
#                           + BG_cons[BG_cons['x_k'] == j]['Q_kj'].sum() 
#                           for j in range(N)])
#
#Q = get_Q(BG_cons)
#
#BG_ints['tau'] = get_taus(Q)
#
#
#u0 = np.zeros(N)
#u0[0] = 1.0
#
#
#print "computing u"
#
#a = 0.0
#b = 2.0
#n = 20000
#
#t = np.linspace(a, b, n)
#u, du = get_us(Q, a, b, n)
#
#
#### CHECK THAT TAU CALCULATION AGREES WITH PDF CALCULATION
##tau_from_Q = np.dot(du,t)/np.sum(du, axis=1)
##
##print BG_ints['tau'] - tau_from_Q
##print BG_ints['tau']/tau_from_Q
##
##
#
#
#
###print "Minimum hitting time for",poly_name
###print "\tbeta =\t", betas[(taus).argmin()]
###print "\ttime =\t", (taus).min()
###print ''
###print "Minimum assembly to disassembly ratio for",poly_name
###print "\tbeta =\t", betas[(taus*nus).argmin()]
###print "\tratio =\t", (taus*nus).min()
##
##print "Minimum hitting time for",poly_name
##print "\talpha =\t", alphas[(taus).argmin()]
##print "\ttime =\t", (taus).min()
##print ''
##print "Minimum assembly to disassembly ratio for",poly_name
##print "\talpha =\t", alphas[(taus*nus).argmin()]
##print "\tratio =\t", (taus*nus).min()
