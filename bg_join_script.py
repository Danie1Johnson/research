import numpy as np

import polyhedra as poly
import manifold_reflected_brownian_motion as mrbm
import bg_run_script as bgrs
import bga_4_0 as bga
import statistics as sts

bga = reload(bga)



def load_ints(poly_name, run_str, output_str):
    """
    Load ensemble of intermediates from files and return a dictionary {int_num: mrbm_obj}
    """
    int_dict = {}
    ints, ids, paths, shell_int, shell_paths, edges, shell_edge = bga.get_bg_ss(poly_name)
    for int_num in range(1, len(ints)):
        filename = bgrs.output_file_name(poly_name, int_num, run_str, output_str)
        int_dict[int_num] = mrbm.MRBM.load(filename)
    return int_dict

def find_sim_edge_relations(int_dict, edges):
    """
    For all of the calculated simulation statistics, find the corresponding bgss edge and
    return a map for relation.
    """
    
    edge_rate_inds = {tuple(e): [] for e in edges}
    
    for int_num, z in int_dict.iteritems():
        if z.hist.hist.shape[0] == 0 or z.hist.hist.sum() == 0.0:
            continue

        int_faces = ints[int_num]
        verts_new, faces_new, face_inds_new = bga.reindex_vertices(face_inds, V, dual, int_faces, verts)
        stat_faces, stat_verts = sts.get_bg_stat_info(int_faces, adj_list, face_inds, face_inds_new)

        for k in range(len(stat_faces)):
            x_new = np.copy(int_faces)
            #print int_faces, x_new, stat_faces
            assert x_new[stat_faces[k][0]] == 0, "Added face already present"
            x_new[stat_faces[k][0]] = 1

            new_int_num = bga.find_int_num(x_new, ints, Rs)
            assert new_int_num != None, "Connection not found"
            
            edge_rate_inds[(int_num, new_int_num)].append(k)

    return edge_rate_inds

def transision_energies(Es, edges):
    """
    Return E_j the value of Es at each of the first int in edges and E_k the second 
    """
    E_j = np.zeros(len(edges))    
    E_k = np.zeros(len(edges))

    for i, edge in enumerate(edges):
        E_j[i] = Es[edge[0]]
        E_k[i] = Es[edge[1]]
    
    return E_j, E_k

def get_rates(epsilon, beta, int_dict, edges, edge_rate_inds, E_j, E_k, Ss, Ts):
    """
    Compute forward and then backward transition rates for bgss.
    """
    forward_rates = np.zeros(len(edges))
    backward_rates = np.zeros(len(edges))

    
    for k, e in enumerate(edges):
        if len(edge_rate_inds[tuple(e)]) == 0:
            forward_rates[k] = 1.0
        else:
            forward_rates[k] = get_rate(epsilon, e[0], edge_rate_inds[tuple(e)], int_dict)
            
    forward_rates *= Ss
    E_jk = E_j - (1.0/beta)*np.log(forward_rates/Ss)
    backward_rates = Ts*np.exp(-beta*(E_jk - E_k))

    return forward_rates, backward_rates, E_jk

def get_rate(epsilon, int_num, stat_inds, int_dict):
    """
    Compute and return the 
    """
    int_rates = int_dict[int_num].hist.rates(np.pi/3.0, epsilon)
    return int_rates[stat_inds].mean()

## for list: invalid value in arccos always in int 7 (or 8)


poly_name = "octahedron"

output_str = "final"

run_str = "test3"

if __name__ == "__main__":

    verts, face_inds, cents = getattr(poly, poly_name)()
    V, E, F, S, species, f_types, adj_list, dual = bga.get_poly(poly_name)
    ints, ids, paths, shell_int, shell_paths, edges, shell_edge = bga.get_bg_ss(poly_name)
    Rs = bga.generate_rotations(adj_list)

    num_ints = len(ints)

    Ss, Ts = bga.get_degeneracies(ints, edges, adj_list, Rs=Rs)
    Es = -bga.get_closed_edges(ints, adj_list)
    E_j, E_k = transision_energies(Es, edges)

    int_dict = load_ints(poly_name, run_str, output_str)
    
    edge_rate_inds = find_sim_edge_relations(int_dict, edges)
    
    epsilon = 0.1
    beta = 1.0

    forward_rates, backward_rates, E_jk = get_rates(epsilon, 
                                              beta, 
                                              int_dict, 
                                              edges, 
                                              edge_rate_inds, 
                                              E_j, 
                                              E_k, 
                                              Ss, 
                                              Ts)



