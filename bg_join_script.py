import numpy as np

import polyhedra as poly
import manifold_reflected_brownian_motion as mrbm
import bg_run_script as bgrs
import bga_4_0 as bga


## for list: invalid value in arccos always in int 7 (or 8)


poly_name = "octahedron"

output_str = "final"

run_str = "test3"

verts, face_inds, cents = getattr(poly, poly_name)()
V, E, F, S, species, f_types, adj_list, dual = bga.get_poly(poly_name)
ints, ids, paths, shell_int, shell_paths, edges, shell_edge = bga.get_bg_ss(poly_name)

num_ints = len(ints)

if __name__ == "__main__":

    for int_num in range(1, num_ints):
        filename = bgrs.output_file_name(poly_name, int_num, run_str, output_str)
        z = mrbm.MRBM.load(filename)
        print z.hist.hist.sum()

