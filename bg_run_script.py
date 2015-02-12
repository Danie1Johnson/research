import datetime
import numpy as np
import platform
import multiprocessing as mp

import bga_4_0 as bga
import manifold_reflected_brownian_motion as mrbm
import triangle_geometry as ti
import matplotlib.pyplot as plt
import polyhedra as poly
bga = reload(bga)
mrbm = reload(mrbm)


def print_sysinfo():

    print '\nPython version  :', platform.python_version()
    print 'compiler        :', platform.python_compiler()

    print '\nsystem     :', platform.system()
    print 'release    :', platform.release()
    print 'machine    :', platform.machine()
    print 'processor  :', platform.processor()
    print 'CPU count  :', mp.cpu_count()
    print 'interpreter:', platform.architecture()[0]
    #print '\n\n'
    print '\n'

def get_time_str():
    dts = str(datetime.datetime.now())
    date_time_str = dts[:10] + "-" + dts[11:13] + "-" + dts[14:16]
    


def run_ints(int_list, total_samples, archive_rate, output_rate, run_str, processor_num):
    """
    Get list of ints to run. Create MRBM instance for each int.
    Run each int for output_rate samples and dump data. Repeat until 
    each int has N samples. Output to a non-overwritten file every
    archive_rate_iterations. 
    """
    
    ### OPTIMIZE PARAMETERS
    int_dict = {}
    for int_num in int_list:
        print "Loading", poly_name, "intermediate", str(int_num) + "..."
        x0, links, lengths, faces = bga.load_bg_int(poly_name, int_num)
        
        # Construct kwargs
        standard_kwarg = {'poly_name': poly_name, 
                          'int_num': int_num}
       
        manifold_kwargs = {'poly_name': poly_name, 
                           'int_num': int_num, 
                           'fixed_com': True,
                           'fixed_rotation': True}
        unary_boundary_kwargs = standard_kwarg
        binary_boundary_kwargs = standard_kwarg 
        stat_kwargs = standard_kwarg
        kwargs = {'manifold_name': manifold_name,  
                  'unary_boundary_name': unary_boundary_name,
                  'binary_boundary_name': binary_boundary_name,
                  'stat_name': stat_name,
                  'manifold_kwargs': manifold_kwargs,
                  'unary_boundary_kwargs': unary_boundary_kwargs,
                  'binary_boundary_kwargs': binary_boundary_kwargs,
                  'stat_kwargs': stat_kwargs,
                  'record_hist': True, 
                  'hist_min': hist_min, 
                  'hist_max': hist_max, 
                  'hist_bins': hist_bins,
                  'err_tol': err_tol}
        
        # Initialize process 
        z = mrbm.MRBM(x0, h0, **kwargs)
        z.optimize_parameters()
        z.reset()
        int_dict[int_num] = z

    ### RUN SIMULATION
    num_samples = 0
    while num_samples < total_samples:
        if total_samples - num_samples < output_rate:
            num_new_samples = total_samples - num_samples
        else:
            num_new_samples = output_rate
        num_samples += num_new_samples
        
        for int_num in int_list:
            int_dict[int_num].sample(N=num_new_samples)
            int_dict.dump(output_file_name(poly_name, int_num, run_str, 'curr'))
            if num_samples == total_samples:
                int_dict.dump(output_file_name(poly_name, int_num, run_str, 'final'))
                print "Intermediate", int_num, "completed." 
            elif num_samples/archive_rate  - (num_samples - num_new_samples)/archive_rate > 0:
                int_dict.dump(output_file_name(poly_name, int_num, run_str, str(num_samples)))
                print "Intermediate", int_num, "archive written with", num_samples, "samples." 
        print "Processor", processor_num, 
        print "completed", num_samples, "samples (",
        print 100*round(float(total_samples)/num_samples, 3), "% )."
    print "Processor", processor_num, "finished."
        

def output_file_name(poly_name, int_num, run_str, output_str):
    filename =  "./results/" + polyname + "/" + poly_name + "_" 
    filename += run_str + "_" + output_str + "_" + str(int_num) + ".csv"
    return filename

### SIMULATION PARAMETERS
manifold_name = 'building_game'

#poly_name = 'tetrahedron'
poly_name = 'octahedron'
#poly_name = 'icosahedron'

unary_boundary_name = 'self_intersection'
binary_boundary_name = 'dihedrals'

stat_name = 'bg_attachment'

verts, face_inds, cents = getattr(poly, poly_name)()
V, E, F, S, species, f_types, adj_list, dual = bga.get_poly(poly_name)
ints, ids, paths, shell_int, shell_paths, edges, shell_edge = bga.get_bg_ss(poly_name)
