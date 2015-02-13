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
    print
    #print '\n\n'
    #print '\n'

def get_time_str():
    dts = str(datetime.datetime.now())
    date_time_str = dts[:10] + "-" + dts[11:13] + "-" + dts[14:16]
    return date_time_str


def run_ints(int_list, total_samples, archive_rate, output_rate, run_str, processor_num):
    """
    Get list of ints to run. Create MRBM instance for each int.
    Run each int for output_rate samples and dump data. Repeat until 
    each int has N samples. Output to a non-overwritten file every
    archive_rate_iterations. 
    """
    hist_min = 0.0
    hist_max = np.pi
    hist_bins = 1000

    err_tol = 10**-12
    h0 = 0.05


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
        #z.optimize_parameters()
        #z.reset()
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
            int_dict[int_num].dump(output_file_name(poly_name, int_num, run_str, 'curr'))
            if num_samples == total_samples:
                int_dict[int_num].dump(output_file_name(poly_name, int_num, run_str, 'final'))
                print "Intermediate", int_num, "completed." 
            elif num_samples/archive_rate  - (num_samples - num_new_samples)/archive_rate > 0:
                int_dict[int_num].dump(output_file_name(poly_name, int_num, run_str, str(num_samples)))
                print "Intermediate", int_num, "archive written with", num_samples, "samples." 
        print "Processor", processor_num, 
        print "completed", num_samples, "samples (",
        print 100*round(num_samples/float(total_samples), 3), "% )."
    print "Processor", processor_num, "finished."
        

def output_file_name(poly_name, int_num, run_str, output_str):
    filename =  "./results/" + poly_name + "/" + poly_name + "_" 
    filename += run_str + "_" + output_str + "_" + str(int_num) + ".pkl"
    return filename


def group_ints(num_groups, int_nums, weights=None):
    """
    Split the ints 0 to len(rates) up into num_groups lists such that each list has 
    approximately the same sum of rates.
    """
    if weights == None:
        weights = np.ones_like(int_nums)

    target = sum(weights)/num_groups*np.arange(num_groups)
    cum_sum = np.cumsum(weights)
    partition = [sum(cum_sum < target[k]) for k in range(num_groups)]
    partition.append(len(weights))
    partition = np.array(partition)
    int_groups = [int_nums[partition[k]:partition[k+1]] for k in range(num_groups)]
    return int_groups
    
    
def get_sample_rate(int_num, h, sample_time):
    return 1.0


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

num_processes = 4

N = 10**2
archive_rate = 50
output_rate = 10

h = 0.5

sample_time = 1.0

run_str = "test2"

if __name__ == "__main__":
    print_sysinfo()
    print get_time_str()

    # Get sample rates, find trivial intermediates
    sample_rates = np.zeros(len(ints))

    for int_num in range(1,len(ints)):
        sample_rates[int_num] = get_sample_rate(int_num, h, sample_time)

    # Split up ints
    int_groups = group_ints(num_processes, np.arange(1,len(ints)), weights=sample_rates)
        
    processes = [mp.Process(target=run_ints, 
                            args=(int_groups[k], 
                                  N, 
                                  archive_rate, 
                                  output_rate, 
                                  run_str, 
                                  k)) for k in range(num_processes)]

    for p in processes:
        p.start()
    for p in processes:
        p.join()
    
