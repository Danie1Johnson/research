import numpy as np
import numpy.linalg

import bga_4_0 as bga
import polyhedra as poly

def get_stat(stat_name, kwargs={}):
    try:
        s_fun = globals()[stat_name](**kwargs)
        return s_fun
    except KeyError, NameError:
        print "ERROR:", stat_name, "not found."
        raise


def bg_attachment(poly_name=None, int_num=None): 
    try:
        # For Building Game intermediates. Denoted 'polyname'.
        assert poly_name != None, "ERROR, no polyhedron specified."
        assert int_num != None, "ERROR, no intermediate specified."
        verts, face_inds, cents = getattr(poly, poly_name)()
        q0, links, lengths, faces = bga.load_bg_int(poly_name, int_num)
        V, E, F, S, species, f_types, adj_list, dual = bga.get_poly(poly_name)
        ints, ids, paths, shell_int, shell_paths, edges, shell_edge = bga.get_bg_ss(poly_name)
        try:
            int_faces = ints[int_num]
            verts_new, faces_new, face_inds_new = bga.reindex_vertices(face_inds, V, dual, int_faces, verts)
            #verts_new, faces_new, face_inds_new = bga.reindex_vertices(face_inds, V, dual, int_faces)
        except IndexError:
            raise Exception("ERROR: " + poly_name + " does not have an intermediate " + int_num)
    except (ValueError, IndexError):
        raise
    
    stat_faces, stat_verts = get_bg_stat_info(int_faces, adj_list, face_inds, face_inds_new)
    
    return lambda x: bg_attachment_fun(x, stat_verts)


def bg_attachment_fun(x, stat_verts):
    """
    Output ndarray with angle for each triplet of vertices in stat_verts
    """

    stat = np.zeros(len(stat_verts))
    for k, vert_trio in enumerate(stat_verts):
        stat[k] = angle_between_edges(x, vert_trio[1], vert_trio[0], vert_trio[2])
    return stat

def get_bg_stat_info(int_faces, adj_list, face_inds, face_inds_new):
    """
    Out put list of faces and list of verts for each stat.
    """
    stat_faces = []
    stat_verts = []

    for k in range(len(int_faces)):
        # Check if face already exists.
        if int_faces[k] != 0:
            continue
        else:
            # See if there are any adjacent faces.
            for j in range(len(adj_list[k])):
                if int_faces[adj_list[k][j]] != 0 and int_faces[adj_list[k][j-1]] != 0:
                    stat_faces.append([k, adj_list[k][j], adj_list[k][j-1]])
                    # Find relevant verticies
                    stat_verts_new = find_vertex_ind(k, 
                                                     adj_list[k][j], 
                                                     adj_list[k][j-1], 
                                                     face_inds, 
                                                     face_inds_new)
                    
                    #remaining_verts = set(face_inds_new[k])
                    #remaining_verts.remove(vert_0)
                    #remaining_verts = list(remaining_verts)
                    #stat_verts_new = [vert_0]
                    #print stat_verts_new, vert_0, remaining_verts, k, j
                    if stat_verts_new != None:
                        stat_faces.append([k, adj_list[k][j], adj_list[k][j-1]])
                        stat_verts.append(stat_verts_new)
                    #assert len(stat_verts_new) == 3, "ERROR: stat_verts incorectly computed"
                    
    return stat_faces, stat_verts

def find_vertex_ind(face, adj_face_1, adj_face_2, face_inds, face_inds_new):
    """
    Take a face and two adjacent faces that share a vertex on face. 
    Find the **new** vertex numbers corresponding to the vertices of two edges 
    (one on each adjacent face) that face meet at the shared vertex. If none, return None.
    """
    ### Get vertex index in case of completed polyhedron 
    original_ind = None
    #for vert in range(V):
    #    if vert in face_inds[face] and vert in race_inds[adj_face_1] and vert in face_inds[adj_face_2]:
    for vert in face_inds[face]:
        if vert in face_inds[adj_face_1] and vert in face_inds[adj_face_2]:
            original_ind = vert
            original_adj_ind_1 = None
            original_adj_ind_2 = None
            for vert_1 in face_inds[face]:
                if vert_1 != vert and vert_1 in face_inds[adj_face_1]:
                    original_adj_ind_1 = vert_1
                    break
            for vert_2 in face_inds[face]:
                if vert_2 != vert and vert_2 in face_inds[adj_face_2]:
                    original_adj_ind_2 = vert_2
                    break
            break
    assert original_ind != None, "Error: common vertex not found"
    assert original_adj_ind_1 != None, "Error: common vertex not found"
    assert original_adj_ind_2 != None, "Error: common vertex not found"

    new_ind_1 = face_inds_new[adj_face_1][face_inds[adj_face_1].index(original_ind)]
    new_ind_2 = face_inds_new[adj_face_2][face_inds[adj_face_2].index(original_ind)]
 
    new_adj_ind_1 = face_inds_new[adj_face_1][face_inds[adj_face_1].index(original_adj_ind_1)]
    new_adj_ind_2 = face_inds_new[adj_face_2][face_inds[adj_face_2].index(original_adj_ind_2)]

    if new_ind_1 == new_ind_2:
        return [new_ind_1, new_adj_ind_1, new_adj_ind_2] 
    else:
        return None


    #print face, adj_face_1, adj_face_2, face_inds, face_inds_new
    #print original_ind
    #### Map index to indexing system for the current intermediate.
    #return face_inds_new[face][face_inds[face].index(original_ind)]
    

def test_1():
    """
    Test for uniformity of dihedral angle in two linked triangles.
    """
    #return lambda x: np.array([dihedral_angle(x, 0, 2, 3, 1)]) 
    return lambda x: np.array([signed_dihedral_angle(x, 0, 2, 3, 1)]) 

def test_2():
    """
    Test for uniformity of dihedral angle in two linked triangles.
    """
    return lambda x: np.array([signed_dihedral_angle(x, 0, 3, 4, 1), 
                               signed_dihedral_angle(x, 3, 0, 4, 2)]) 
        
def angle_between_edges(x, v0, v1, v2):
    """
    Return the angle between edge (v0 v1) and edge (v2 v1). 
    """

    ans = np.dot(x[3*v0:3*v0+3] - x[3*v1:3*v1+3],
                 x[3*v2:3*v2+3] - x[3*v1:3*v1+3])
    ans /= numpy.linalg.norm(x[3*v0:3*v0+3] -x[3*v1:3*v1+3])
    ans /= numpy.linalg.norm(x[3*v2:3*v2+3] -x[3*v1:3*v1+3])
    return np.arccos(ans)

#def dihedral_angle(x, v0, v1a, v1b, v2):
#    """
#    Return the angle between v0, the v1a--v1b midpoint and v2 
#    """
#
#    #ans = np.dot(x[3*v0:3*v0+3] - 0.5*(x[3*v1a:3*v1a+3] + x[3*v1b:3*v1b+3]),
#    #             x[3*v2:3*v2+3] - 0.5*(x[3*v1a:3*v1a+3] + x[3*v1b:3*v1b+3]))
#    #ans /= numpy.linalg.norm(x[3*v0:3*v0+3] - 0.5*(x[3*v1a:3*v1a+3] + x[3*v1b:3*v1b+3]))
#    #ans /= numpy.linalg.norm(x[3*v2:3*v2+3] - 0.5*(x[3*v1a:3*v1a+3] + x[3*v1b:3*v1b+3]))
#    #return np.arccos(ans)


#def find_dihedrals(f1, f2, x, faces):
#    """
#    Take the index of two faces of a BG intermediate and compute the dihedral angle between them
#    """
#    n1 = triangle_normal(f1, x, faces)
#    n2 = triangle_normal(f2, x, faces)
#    return np.arccos(np.dot(n1, n2))
#
#def triangle_normal(f, x, faces):
#    """
#    Find normal vector to triangle. 
#    """
#    v1 = x[3*faces[f][0]:3*faces[f][0]+3] - x[3*faces[f][1]:3*faces[f][1]+3] 
#    v2 = x[3*faces[f][2]:3*faces[f][2]+3] - x[3*faces[f][1]:3*faces[f][1]+3] 
#    n = np.cross(v1, v2)
#    return n/numpy.linalg.norm(n)


def triangle_normal(x, va, vb, vc):
    """
    Find the normal vector at va
    """ 
    return np.cross(x[3*vb:3*vb+3] - x[3*va:3*va+3], x[3*vc:3*vc+3] - x[3*va:3*va+3]) 

def signed_dihedral_angle2(x, v0, v1a, v1b, v2):
    """
    Return the angle between v0, the v1a--v1b midpoint and v2 
    """
    # Triangle normals at v1a
    n1 = triangle_normal(x, v1a, v1b, v0)
    n2 = triangle_normal(x, v1a, v1b, v2)

    ans = np.arccos(np.dot(n1, n2)/(numpy.linalg.norm(n1)*numpy.linalg.norm(n2)))
    
    ###
    ans = np.pi - ans

    if np.dot(n1, x[3*v2:3*v2+3]- x[3*v1a:3*v1a+3]) > 0:
        return ans
    else:
        return 2.0*np.pi - ans

def signed_dihedral_angle(x, v0, v1a, v1b, v2):
    """
    Return the angle between v0, the v1a--v1b midpoint and v2 
    """

    ans = np.dot(x[3*v0:3*v0+3] - 0.5*(x[3*v1a:3*v1a+3] + x[3*v1b:3*v1b+3]),
                 x[3*v2:3*v2+3] - 0.5*(x[3*v1a:3*v1a+3] + x[3*v1b:3*v1b+3]))
    ans /= numpy.linalg.norm(x[3*v0:3*v0+3] -0.5*(x[3*v1a:3*v1a+3] + x[3*v1b:3*v1b+3]))
    ans /= numpy.linalg.norm(x[3*v2:3*v2+3] -0.5*(x[3*v1a:3*v1a+3] + x[3*v1b:3*v1b+3]))

    n = np.cross(x[3*v1a:3*v1a+3] - x[3*v0:3*v0+3], x[3*v1b:3*v1b+3] - x[3*v0:3*v0+3]) 
    theta = np.arccos(ans)
    if np.dot(x[3*v2:3*v2+3] - x[3*v0:3*v0+3], n) > 0.0:
        return theta
    else:
        return 2.0*np.pi - theta
    
