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

import mpmath as mp
import matplotlib.pyplot as plt
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
    

def get_bg_ss(poly_name, get_degens=False):
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
    if get_degens == True:
        degens = []
    for j in range(E):
        line = [int(x) for x in f.readline().split()]
        edges.append(line[0:2])
        shell_edge.append(line[2])
        if get_degens == True:
            degens.append(line[3])

    edges = np.array(edges)
    shell_edge = np.array(shell_edge)


    if get_degens == True:
        return ints, ids, paths, shell_int, shell_paths, edges, shell_edge, degens
    else:
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
    
    

def edge_adj_list(edges, inds=False):
    """
    Take list of edge pairs and populate a forward adjacency list.
    """
    num_ints = 1 + max(max(e) for e in edges)
    
    edge_adj_list = [[] for k in range(num_ints)]
    edge_adj_inds = [[] for k in range(num_ints)]

    for k, e in enumerate(edges):
        edge_adj_list[min(e[0], e[1])].append(max(e[0],e[1]))
        edge_adj_inds[min(e[0], e[1])].append(k)

    if inds == True:
        return edge_adj_list, edge_adj_inds
    else:
        return edge_adj_list

def get_paths(edges, ints, shell_int=None, shellable=False):
    """
    Find number of pathways to each intermediate.
    """

    mp.mp.dps = 100
    
    int_sizes = [sum(np.array(int_j) != 0) for int_j in ints]
    edge_adj = edge_adj_list(edges)

    # Initialize single faced ints to 1 path
    paths = [mp.mpf(int(sum(np.array(int_j) != 0) == 1)) for int_j in ints]
    F = len(ints[0])
    for n in range(1,F):
        for i in range(len(ints)):
            if int_sizes[i] == n:
                for k in edge_adj[i]:
                    #print n, i, k
                    if shellable == True:
                        if shell_int[i] != 0 and shell_int[k] != 0:
                            paths[k] += paths[i]
                    else:
                        paths[k] += paths[i]
                
    return paths


def get_shellings(edges, ints, Ss, shell_int, poly_adj_list):
    """
    Find number of pathways to each intermediate.
    """
    mp.mp.dps = 100

    Rs = Rs = generate_rotations(poly_adj_list)
    int_sizes = [sum(np.array(int_j) != 0) for int_j in ints]
    edge_adj, edge_adj_inds = edge_adj_list(edges, inds=True)
    
    # Initialize single faced ints 
    #shells = [int(sum(np.array(int_j) != 0) == 1) for int_j in ints]
    shells = [mp.mpf(0.0) for x in ints]
    for m in range(len(ints)):
        if int_sizes[m] == 1:
            shells[m] += len(Rs)/get_r(Rs, ints[m])
    
    F = len(ints[0])
    for n in range(1,F):
        for i in range(len(ints)):
            if int_sizes[i] == n:
                for ind, k in enumerate(edge_adj[i]):
                    if shell_int[i] != 0 and shell_int[k] != 0:
                        #print 'hi'
                        shells[k] += shells[i]*Ss[edge_adj_inds[i][ind]]
                
    return shells

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

def get_closed_edges(ints, adj_list):
    """
    For each intermediate in ints, return the number of closed edges it has.
    """
    Es = []
    
    for inter in ints:
        E = 0
        for k in range(len(inter)):
            if inter[k] == 1:
               for j in adj_list[k]:
                   if inter[j] == 1 and j < k:
                       E += 1
        Es.append(E)

    return np.array(Es)


def get_degeneracies(ints, edges, adj_list, Rs=None):
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


#def generate_rotations(adj_list):
#    """
#    Generate a list of index permutations corresponding to the 
#    rotation group of adj_list's polyhedron.
#    """
#    rotations = []
#    
#    F = len(adj_list)
#    
#    for j in range(F):
#        for k in range(len(adj_list[j])):
#            ind = rotate_ind(j,k,adj_list)
#            if ind == None:
#                continue
#            rotations.append(ind)
#
#    return rotations
#
#def rotate_ind(piv,r,adj_list):
#    """
#    Get index permutation coresponding to moving face 0 to piv 
#    and rotating that face in adj_list r times.
#    """
#    f = 0
#    F = len(adj_list)
#    
#    prev = []
#    nxt = []
#    ind = [None for x in range(F)]
#    
#    if piv >= F:
#        raise Exception("bad piv!")
#    
#    # For reindexing, recenter at piv
#    ind[0] = piv
#    f += 1
#
#    # Number of faces adjacent to piv (must be the same for [0] and [piv])
#    N = len(adj_list[0]) 
#    
#    # If base faces different, return error
#    if len(adj_list[0]) != len(adj_list[piv]):
#        return None
#    
#    # Fill in indices for faces adjacent to piv
#    for j in range(N):
#        ind[adj_list[0][j]] = adj_list[piv][(r+j-1+N)%N]
#        prev.append(adj_list[0][j])
#        f += 1
#
#    while f < F:
#        if len(prev) == 0:
#            print "ERROR, inf-loop"
#            break
#        
#        for k in range(len(prev)):
#            Nk = len(adj_list[prev[k]])
#
#            for h in range(Nk):
#                # Find already indexed face
#                i = adj_list[prev[k]][h]
#                if ind[i] != -1:
#                    break
#                
#            for q in range(Nk):
#                if adj_list[ind[prev[k]]][q] == ind[i]:
#                    break
#
#            for g in range(Nk):
#                map_from = adj_list[prev[k]][(h+g)%Nk]
#                map_to = adj_list[ind[prev[k]]][(q+g)%Nk]
#
#                # Check if already indexed 
#                if ind[map_from] != -1:
#                    if ind[map_from] != map_to:
#                        return None
#
#                else:
#                    # If not, add index
#                    ind[map_from] = map_to
#                    nxt.append(map_from)
#                    f += 1
#              
#        prev = []
#        prev = copy.copy(nxt)
#        nxt = []
#    return ind
#

def generate_rotations(adj_list):
    """
    Generate a list of index permutations corresponding to the 
    rotation group of adj_list's polyhedron.
    """
    rotations = []
    
    F = len(adj_list)
    
    for j in range(F):
        for k in range(len(adj_list[j])):
            try:
                ind = rotate_ind(j,k,adj_list)
                if ind[0] == -2:
                    continue
            except:
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
                    #print 'a'
                    break
                
            for q in range(Nk):
                if adj_list[ind[prev[k]]][q] == ind[i]:
                    #print 'b'
                    break

            for g in range(Nk):
                map_from = adj_list[prev[k]][(h+g)%Nk]
                #print k, g, Nk, ind[prev[k]], len(adj_list), len(adj_list[ind[prev[k]]]), (q+g)%Nk 
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
    if sum(x > 0 for x in x1) + 1 != sum(y > 0 for y in x2):
        print 'ERROR: non-coresponding x1 and x2, wrong number of faces.' 
        print x1, x2

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
    #if numpy.linalg.norm(x1 + x2) == (4*sum(x1)+1)**0.5:
    #    return True
    if sum((x1[k] != 0) != (x2[k] != 0) for k in range(len(x1))) == 1:
        return True
    return False


def same_int(x1,x2,Rs):
    """
    Check if x1 is a polyhedral rotation of x2.
    """
    for R in Rs:
        #if numpy.linalg.norm(np.array([x1[R[k]] for k in range(len(R))]) - np.array(x2)) == 0:
        #    return True
        if sum((x1[R[k]] != 0) != (x2[k] != 0) for k in range(len(x2))) == 0:
            return True
    return False

def find_int_num(x, ints, Rs):
    """
    For binary vector x of faces present determine which bg int (if any)  x belongs to. 
    """
    num_faces = sum(x >= 1)
    for k, int_k in enumerate(ints):
        if sum(int_k >= 1) != num_faces:
            continue
        else:
            if same_int(x, int_k, Rs) == True:
                return k
    return None



def get_rs(ints,adj_list,Rs=None):
    """
    Compute the order of the rotation group for each intermediate.
    """
    if Rs == None:
        Rs = generate_rotations(adj_list)
    rs = []

    for inter in ints:
        #r = 0
        #for R in Rs:
        #    if numpy.linalg.norm(inter - np.array([inter[R[k]] for k in range(len(inter))])) == 0:
        #        r += 1
        #rs.append(r)
        rs.append(get_r(Rs, inter))

    return np.array(rs)

def get_r(Rs, inter):
    r = 0
    for R in Rs:
        if numpy.linalg.norm(inter - np.array([inter[R[k]] for k in range(len(inter))])) == 0:
            r += 1
    return r


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

def get_dist(Q, T, num_times):
    """
    Return 2d array with analytic solution for the dist on each intermediate starting from int 0.
    """
    N = Q.shape[0]

    u0 = np.zeros(N)
    u0[0] = 1.0

    # Identity matrix excpt with M_j,j = 0 for if j in A
    Ident_Ac = scipy.sparse.coo_matrix((np.ones(N-1),(np.arange(N-1),np.arange(N-1))), shape=(N,N))

    Ident_Ac = Ident_Ac.tocsc()
    Q = Q.tocsc()
    
    P = scipy.sparse.linalg.expm_multiply(Q, 
                                          u0, 
                                          start=0.0, 
                                          stop=T, 
                                          num=num_times, 
                                          endpoint=True)

    return P

def get_taus(Q, N, A=None):
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
    


def get_connected_faces(dual_v, int_faces):
    """
    Given a list of face indices dual_v of faces meeting at a particular vertex
    and the inclussion/exclussion vector int_faces indicating which faces are in the int,
    return a set of sets where each inner set contains the face indices that are edge connected 
    at the vertex.
    """
    face_groups = set()
    for f in dual_v:
        if int_faces[f] != 0:
            face_groups.add(frozenset([f]))
    for k, f_1 in enumerate(dual_v):
        f_0 = dual_v[k-1]
        if int_faces[f_0] != 0 and int_faces[f_1] != 0:
            fg_0 = None
            fg_1 = None
            for fg in face_groups:
                if f_0 in fg:
                    fg_0 = fg
                if f_1 in fg:
                    fg_1 = fg
            if fg_0 == fg_1:
                continue
            try:
                fg_new = fg_0.union(fg_1)
                face_groups.remove(fg_0)
                face_groups.remove(fg_1)
                face_groups.add(fg_new)
            except KeyError:
                print "ERROR:", fg_0, "or", fg_1, "not found in", face_groups
                raise
    return face_groups

def verify_dual_order(dual):
    ###### NEEDS ADDING?!?!?! ###########
    return dual
    
def load_bg_int(poly_name, int_num):
    """
    Load information about specified polyhedral intermediate. 
    Return information for corresponding linkage.
    """

    try: 
        poly_info = getattr(poly, poly_name)
    except AttributeError:
        raise Exception("ERROR: " + poly_name + " not found in polyhedra.py")


    verts, face_inds, cents = poly_info()
    V, E, F, S, species, f_types, adj_list, dual = get_poly(poly_name)
    ints, ids, paths, shell_int, shell_paths, edges, shell_edge = get_bg_ss(poly_name)
    dual = verify_dual_order(dual)

    try:
        int_faces = ints[int_num]
    except IndexError:
        raise Exception("ERROR: " + poly_name + " does not have an intermediate " + int_num)
        
    # Reindex vertices/
    verts_new, faces, face_inds_new = reindex_vertices(face_inds, V, dual, int_faces, verts)
    x0 = verts_new.flatten()

    #face_inds_new = [[None for f in fi] for fi in face_inds]
    #verts_new = []
    #new_V = 0
    #
    #for v in range(V):
    #    f_groups = get_connected_faces(dual[v], int_faces)
    #    for fg in f_groups:
    #        verts_new.append(verts[v,:])
    #        for f in fg:
    #            face_inds_new[f][face_inds[f].index(v)] = new_V
    #        new_V += 1
    #verts = np.array(verts_new)
    #x0 = verts.flatten()
    #
    ##print 'A', face_inds_new
    #faces = [f for k, f in enumerate(face_inds_new) if int_faces[k] != 0]

    #print 'B',faces
    # Make list of links.
    links = set()
    for face in faces:
        #print face
        for k in range(len(face)):
            #print '\t', k, face[k], face[k-1]
            links.add(frozenset([face[k], face[k-1]]))
    links = [list(link) for link in links]
    lengths = np.array([numpy.linalg.norm(verts_new[link[0],:] - verts_new[link[1],:]) for link in links])
    
    return x0, links, lengths, faces


def reindex_vertices(face_inds, V, dual, int_faces, verts):
    """
    """
    face_inds_new = [[None for f in fi] for fi in face_inds]
    verts_new = []
    new_V = 0

    for v in range(V):
        f_groups = get_connected_faces(dual[v], int_faces)
        for fg in f_groups:
            verts_new.append(verts[v,:])
            for f in fg:
                face_inds_new[f][face_inds[f].index(v)] = new_V
            new_V += 1
    verts_new = np.array(verts_new)
    x0 = verts_new.flatten()

    #print 'A', face_inds_new
    faces = [f for k, f in enumerate(face_inds_new) if int_faces[k] != 0]
    return verts_new, faces, face_inds_new

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
    bg_int.sample(progress_bar=False)
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


def plot_polyhedron(v, f_inds, inter=np.array([])):
    if len(inter) == 0:
        inter = np.ones((len(f_inds),1))
    ax = Axes3D(plt.figure())
    scale = np.abs(v).max()*1.2
    ax.set_xlim(-scale,scale)
    ax.set_ylim(-scale,scale)
    ax.set_zlim(-scale,scale)
    for i in range(len(f_inds)):
        if inter[i] != 0:
            side = []
            for j in range(len(f_inds[i])):
                #side1 = [[0.,0.,0.],[0.,0.,1.],[0.,1.,1.],[0.,1.,0.]]
                side.append([v[f_inds[i][j],0],v[f_inds[i][j],1],v[f_inds[i][j],2]])

            tri = Poly3DCollection([side])
            color = colors.rgb2hex(sp.rand(3))
            tri.set_facecolor(color)
            tri.set_edgecolor('k')

            ax.add_collection3d(tri)

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


def first_face(inter):
    for k,x in enumerate(inter):
        if x != 0:
            return k
    return -1


def face_constraints(v,i0,i1,i2,ik,g0,g1,g2,gk,N):
    #N = 3*v.shape[0]
    # Template vertices
    t0 = v[i0,:].T
    t1 = v[i1,:].T
    t2 = v[i2,:].T
    tk = v[ik,:].T
    # Face vertices
    v0 = v[i0,:].T
    v1 = v[i1,:].T
    v2 = v[i2,:].T
    vk = v[ik,:].T
    # Constants
    gamma = 1.0/numpy.linalg.norm(np.cross(v0-v1,v2-v1))
    tau = numpy.linalg.norm(tk-t1)/numpy.linalg.norm(t0-t1)
    theta = math.acos(np.dot(tk-t1,t0-t1)/(numpy.linalg.norm(tk-t1)*numpy.linalg.norm(t0-t1)))
    #theta = -acos(np.dot(tk-t1,t0-t1)/(numpy.linalg.norm(tk-t1)*numpy.linalg.norm(t0-t1)))
    # Normal
    n = np.cross(v0-v1,v2-v1)*gamma

    #print "cross calc", v0-v1,v2-v1,gamma
    # Rotation matrix
    R = np.zeros((3,3))

    R[0,0] = np.cos(theta) + n[0]**2*(1-np.cos(theta))
    R[1,0] = n[1]*n[0]*(1-np.cos(theta)) + n[2]*np.sin(theta)
    R[2,0] = n[2]*n[0]*(1-np.cos(theta)) - n[1]*np.sin(theta)

    R[0,1] = n[0]*n[1]*(1-np.cos(theta)) - n[2]*np.sin(theta)
    R[1,1] = np.cos(theta) + n[1]**2*(1-np.cos(theta))
    R[2,1] = n[2]*n[1]*(1-np.cos(theta)) + n[0]*np.sin(theta)

    R[0,2] = n[0]*n[2]*(1-np.cos(theta)) + n[1]*np.sin(theta)
    R[1,2] = n[1]*n[2]*(1-np.cos(theta)) - n[0]*np.sin(theta)
    R[2,2] = np.cos(theta) + n[2]**2*(1-np.cos(theta))

    #rotate to see if got the right rotation
    #print "v0,v1,v2,vk",v0,v1,v2,vk
    #print "gamma, theta, tau",gamma,theta,tau
    #print "normal",n
    #vk_hat = v1 + np.dot(R.T,v0-v1)*tau
    vk_hat = v1 + np.dot(R,v0-v1)*tau
    #print "vks",vk,vk_hat

    # Normal Jacobian
    K = np.zeros((3,9))

    K[0,0] = gamma*0
    K[1,0] = gamma*(-v2[2]+v1[2])
    K[2,0] = gamma*(+v2[1]-v1[1])

    K[0,1] = gamma*(+v2[2]-v1[2])
    K[1,1] = gamma*0
    K[2,1] = gamma*(-v2[0]+v1[0])

    K[0,2] = gamma*(-v2[1]+v1[1])
    K[1,2] = gamma*(+v2[0]-v1[0])
    K[2,2] = gamma*0

    K[0,3] = gamma*0
    K[1,3] = gamma*(+v2[2]-v0[2])
    K[2,3] = gamma*(-v2[1]+v0[1])

    K[0,4] = gamma*(-v2[2]+v0[2])
    K[1,4] = gamma*0
    K[2,4] = gamma*(+v2[0]-v0[0])

    K[0,5] = gamma*(+v2[1]-v0[1])
    K[1,5] = gamma*(-v2[0]+v0[0])
    K[2,5] = gamma*0

    K[0,6] = gamma*0
    K[1,6] = gamma*(+v0[2]-v1[2])
    K[2,6] = gamma*(-v0[1]+v1[1])

    K[0,7] = gamma*(-v0[2]+v1[2])
    K[1,7] = gamma*0
    K[2,7] = gamma*(+v0[0]-v1[0])

    K[0,8] = gamma*(+v0[1]-v1[1])
    K[1,8] = gamma*(-v0[0]+v1[0])
    K[2,8] = gamma*0

    new_row = np.zeros((3,N))

    new_row[:,3*g0+0] = tau*(np.dot(R_partial(n,K[:,0],theta),v0-v1) + R[:,0])
    new_row[:,3*g0+1] = tau*(np.dot(R_partial(n,K[:,1],theta),v0-v1) + R[:,1])
    new_row[:,3*g0+2] = tau*(np.dot(R_partial(n,K[:,2],theta),v0-v1) + R[:,2])

    new_row[:,3*g1+0] = tau*(np.dot(R_partial(n,K[:,3],theta),v0-v1) - R[:,0]) + np.array([1.0,0,0]).T 
    new_row[:,3*g1+1] = tau*(np.dot(R_partial(n,K[:,4],theta),v0-v1) - R[:,1]) + np.array([0,1.0,0]).T 
    new_row[:,3*g1+2] = tau*(np.dot(R_partial(n,K[:,5],theta),v0-v1) - R[:,2]) + np.array([0,0,1.0]).T 

    new_row[:,3*g2+0] = tau*(np.dot(R_partial(n,K[:,6],theta),v0-v1))
    new_row[:,3*g2+1] = tau*(np.dot(R_partial(n,K[:,7],theta),v0-v1))
    new_row[:,3*g2+2] = tau*(np.dot(R_partial(n,K[:,8],theta),v0-v1))

    new_row[:,3*gk+0] = np.array([-1.0,0,0]).T
    new_row[:,3*gk+1] = np.array([0,-1.0,0]).T
    new_row[:,3*gk+2] = np.array([0,0,-1.0]).T

    return new_row


def hasedge(a,b,f_inds):
    #a,b are faces
    count = 0
    for v in f_inds[a]: 
        if v in f_inds[b]:
            count += 1
    if count == 2:
        return True
    return False

def R_partial(n,dn,theta):
    Rp = np.zeros((3,3))
    
    Rp[0,0] = 2.0*n[0]*dn[0]*(1.0-np.cos(theta))
    Rp[1,0] = (1.0-np.cos(theta))*(n[1]*dn[0] + dn[1]*n[0]) + dn[2]*np.sin(theta)
    Rp[2,0] = (1.0-np.cos(theta))*(n[2]*dn[0] + dn[2]*n[0]) - dn[1]*np.sin(theta)
    
    Rp[0,1] = (1.0-np.cos(theta))*(n[0]*dn[1] + dn[0]*n[1]) - dn[2]*np.sin(theta)
    Rp[1,1] = 2.0*n[1]*dn[1]*(1.0-np.cos(theta))
    Rp[2,1] = (1.0-np.cos(theta))*(n[2]*dn[1] + dn[2]*n[1]) + dn[0]*np.sin(theta)
    
    Rp[0,2] = (1.0-np.cos(theta))*(n[0]*dn[2] + dn[0]*n[2]) + dn[1]*np.sin(theta)
    Rp[1,2] = (1.0-np.cos(theta))*(n[1]*dn[2] + dn[1]*n[2]) - dn[0]*np.sin(theta)
    Rp[2,2] = 2.0*n[2]*dn[2]*(1.0-np.cos(theta))
    
    return Rp

def DoF(inter,v,f_inds):
    g_inds = []
    #g_inv_inds = []
    g_i = 0
    for q in range(len(inter)):
        g_inds.append(range(g_i,g_i+len(f_inds[q])))
        g_i += len(f_inds[q])
    
    if sum(inter) == 0:
        return 0, np.array([])
                       
    # Base Constraints
    N = 3*(g_i)
    #N = 3*v.shape[0]
    base = first_face(inter)
    if base == -1:
        return 0,np.array([])
    J = np.zeros((3*len(f_inds[base]),N))
    for i in range(len(f_inds[base])):
        J[3*i+0,3*g_inds[base][i]+0] = 1.0
        J[3*i+1,3*g_inds[base][i]+1] = 1.0
        J[3*i+2,3*g_inds[base][i]+2] = 1.0
    # Length Constraints
    for j in range(len(f_inds)):
        if inter[j] != 0:
            # Same Vertex Constraints
            for w in range(j):
                if inter[w] != 0 and hasedge(j,w,f_inds): ##should give same results as before. until this line used
                #if inter[w] != 0:
                    for r in range(len(f_inds[j])):
                        if f_inds[j][r] in f_inds[w]:
                            #w_ind = np.where(f_inds[j][r] == f_inds[w])
                            w_ind = f_inds[w].index(f_inds[j][r])
                            new_row = np.zeros((3,N))
                            new_row[0,3*g_inds[j][r]+0] = 1.0
                            new_row[0,3*g_inds[w][w_ind]+0] = -1.0                          
                            new_row[1,3*g_inds[j][r]+1] = 1.0
                            new_row[1,3*g_inds[w][w_ind]+1] = -1.0                          
                            new_row[2,3*g_inds[j][r]+2] = 1.0
                            new_row[2,3*g_inds[w][w_ind]+2] = -1.0                          
                            J = np.vstack((J,new_row))
                            #print 'SVC'
            # Edge Length Constraints
            for k in range(len(f_inds[j])):
                new_row = np.zeros((1,N))
                new_row[0,3*g_inds[j][k-0]+0] = 2.0*(v[f_inds[j][k-0],0] - v[f_inds[j][k-1],0])                       
                new_row[0,3*g_inds[j][k-1]+0] = 2.0*(v[f_inds[j][k-1],0] - v[f_inds[j][k-0],0])                       
                new_row[0,3*g_inds[j][k-0]+1] = 2.0*(v[f_inds[j][k-0],1] - v[f_inds[j][k-1],1])                       
                new_row[0,3*g_inds[j][k-1]+1] = 2.0*(v[f_inds[j][k-1],1] - v[f_inds[j][k-0],1])                       
                new_row[0,3*g_inds[j][k-0]+2] = 2.0*(v[f_inds[j][k-0],2] - v[f_inds[j][k-1],2])                       
                new_row[0,3*g_inds[j][k-1]+2] = 2.0*(v[f_inds[j][k-1],2] - v[f_inds[j][k-0],2])                       
                J = np.vstack((J,new_row))
            # Face Angle Constraints
            new_row = np.zeros((1,N))
            new_row[0,3*g_inds[j][0]+0] = v[f_inds[j][2],0] - v[f_inds[j][1],0]                          
            new_row[0,3*g_inds[j][1]+0] = 2.0*v[f_inds[j][1],0] - v[f_inds[j][0],0] - v[f_inds[j][2],0] 
            new_row[0,3*g_inds[j][2]+0] = v[f_inds[j][0],0] - v[f_inds[j][1],0]                              
            new_row[0,3*g_inds[j][0]+1] = v[f_inds[j][2],1] - v[f_inds[j][1],1]                          
            new_row[0,3*g_inds[j][1]+1] = 2.0*v[f_inds[j][1],1] - v[f_inds[j][0],1] - v[f_inds[j][2],1] 
            new_row[0,3*g_inds[j][2]+1] = v[f_inds[j][0],1] - v[f_inds[j][1],1]                          
            new_row[0,3*g_inds[j][0]+2] = v[f_inds[j][2],2] - v[f_inds[j][1],2]                          
            new_row[0,3*g_inds[j][1]+2] = 2.0*v[f_inds[j][1],2] - v[f_inds[j][0],2] - v[f_inds[j][2],2] 
            new_row[0,3*g_inds[j][2]+2] = v[f_inds[j][0],2] - v[f_inds[j][1],2]                          
            J = np.vstack((J,new_row))
            for k in range(3,len(f_inds[j])):
                J = np.vstack((J,
                               face_constraints(v,
                                                f_inds[j][0],
                                                f_inds[j][1],
                                                f_inds[j][2],
                                                f_inds[j][k],
                                                g_inds[j][0],
                                                g_inds[j][1],
                                                g_inds[j][2],
                                                g_inds[j][k],
                                                N)))
                
    free_vars = 0
    for k in range(v.shape[0]):
        af = 0
        for i in range(len(f_inds)):
            if k in f_inds[i] and inter[i] != 0:
                af += 1
                break
        if af == 0:
            free_vars += 3

    rank1 = numpy.linalg.matrix_rank(J)
    #u, s, v = np.linalg.svd(J)
    #print u
    #rank2 = np.sum(s > 1e-4)
    #rank2 = -1
    #if rank1 != rank2:
    #     print "oh dear. the ranks is fucked."
    #print 'J size',J.shape,'free variables',free_vars,'rank',rank1,rank2
    #print J
    rank = rank1
    #return N - rank - free_vars, J
    fv = sum(np.array([len(g_inds[i]) for i in range(len(inter)) if inter[i] == 0]))
    #print "returns",N,rank,fv
    return N - rank - 3*fv, J





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
