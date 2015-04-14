import datetime
import numpy as np
import platform
import multiprocessing as mp
import time

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


def run_ints(poly_name, 
             int_list, 
             total_samples, 
             archive_rate, 
             output_rate, 
             run_str, 
             processor_num, 
             h, 
             err_tol,
             skip_sampling):
    """
    Get list of ints to run. Create MRBM instance for each int.
    Run each int for output_rate samples and dump data. Repeat until 
    each int has N samples. Output to a non-overwritten file every
    archive_rate_iterations. 
    """
    ### OPTIMIZE PARAMETERS
    bg_kwargs = get_bg_kwargs(poly_name, int_list, err_tol)
    int_dict = {}
    for int_num in int_list:
        print "Loading", poly_name, "intermediate", str(int_num)
        x0, links, lengths, faces = bga.load_bg_int(poly_name, int_num)
        
        ## Construct kwargs
        #standard_kwarg = {'poly_name': poly_name, 
        #                  'int_num': int_num}
        #
        #manifold_kwargs = {'poly_name': poly_name, 
        #                   'int_num': int_num, 
        #                   'fixed_com': True,
        #                   'fixed_rotation': True}
        #unary_boundary_kwargs = standard_kwarg
        #binary_boundary_kwargs = standard_kwarg 
        #stat_kwargs = standard_kwarg
        #kwargs = {'manifold_name': manifold_name,  
        #          'unary_boundary_name': unary_boundary_name,
        #          'binary_boundary_name': binary_boundary_name,
        #          'stat_name': stat_name,
        #          'manifold_kwargs': manifold_kwargs,
        #          'unary_boundary_kwargs': unary_boundary_kwargs,
        #          'binary_boundary_kwargs': binary_boundary_kwargs,
        #          'stat_kwargs': stat_kwargs,
        #          'record_hist': True, 
        #          'hist_min': hist_min, 
        #          'hist_max': hist_max, 
        #          'hist_bins': hist_bins,
        #          'err_tol': err_tol}
        
        
        # Initialize process 
        z = mrbm.MRBM(x0, h, **bg_kwargs[int_num])
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
            if skip_sampling[int_num] == False:
                int_dict[int_num].sample(N=num_new_samples, record_trace=False, record_stats=False)
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
        

def get_bg_kwargs(poly_name, 
                  int_list, 
                  err_tol, 
                  newton_max_itr=100,
                  fixed_com=True, 
                  fixed_rotation=True, 
                  record_hist=True, 
                  hist_min=0.0, 
                  hist_max=np.pi, 
                  hist_bins=999,
                  manifold_name='building_game',
                  unary_boundary_name='self_intersection',
                  binary_boundary_name='dihedrals',
                  stat_name='bg_attachment'):
    bg_kwargs = {}
    for int_num in int_list:
        x0, links, lengths, faces = bga.load_bg_int(poly_name, int_num)
        
        # Construct kwargs
        standard_kwarg = {'poly_name': poly_name, 
                          'int_num': int_num}
       
        manifold_kwargs = {'poly_name': poly_name, 
                           'int_num': int_num, 
                           'fixed_com': fixed_com,
                           'fixed_rotation': fixed_rotation}
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
                  'record_hist': record_hist, 
                  'hist_min': hist_min, 
                  'hist_max': hist_max, 
                  'hist_bins': hist_bins,
                  'err_tol': err_tol}
        bg_kwargs[int_num] = kwargs
    return bg_kwargs
        

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
    
def get_sample_rates(poly_name, int_list, h, err_tol, samples):
    """
    Sample each int of poly_name in int_list for specified number of samples and 
    return an np array with the time it took.
    """
    print "Drawing samples to estimate sampling rates."
    rates = np.zeros(len(int_list))
    skip_sampling = {int_num: True for int_num in int_list}
    int_dict = {}
    bg_kwargs = get_bg_kwargs(poly_name, int_list, err_tol)
        
    for k, int_num in enumerate(int_list):
        x0, links, lengths, faces = bga.load_bg_int(poly_name, int_num)
        z = mrbm.MRBM(x0, h, **bg_kwargs[int_num])

        if z.num_stats == 0:
            continue
        if z.m <= 0:
            continue
        if z.m - z.manifold_mod_directions(z.x).shape[0] <= 0:
            continue

        skip_sampling[int_num] = False
        before = time.clock()
        z.sample(N=samples, record_trace=False, record_stats=False)
        after = time.clock()
        rates[k] = after - before
        print int_num, z.m - z.manifold_mod_directions(z.x).shape[0], after - before
    print '\nEstimated serial time for 1,000,000 samples:', 10**6*sum(rates)/samples/3600.0, "hours\n"
    return rates, skip_sampling

def optimize_sampling(poly_name, num_ints, err_tol):
    hs = np.arange(0.030,0.065,0.005)
    samples_per_trial = 10**3
    int_list = range(1,num_ints)
    bg_kwargs = get_bg_kwargs(poly_name, int_list, err_tol, newton_max_itr=10000)
    times = np.zeros_like(hs)
    
    int_dict = {}
    for int_num in range(1,num_ints):
       x0, links, lengths, faces = bga.load_bg_int(poly_name, int_num)
       int_dict[int_num] = mrbm.MRBM(x0, hs[0], **bg_kwargs[int_num])
        
    for k, h in enumerate(hs):
        print h
        for int_num in range(1,num_ints):
            if int_dict[int_num].num_stats == 0:
                continue
            if int_dict[int_num].m <= 0:
                continue
            if (int_dict[int_num].m - 
                int_dict[int_num].manifold_mod_directions(int_dict[int_num].x).shape[0]) <= 0:
                continue
            int_dict[int_num].x = int_dict[int_num].x0
            int_dict[int_num].h = h
            before = time.clock()
            int_dict[int_num].sample(N=samples_per_trial, record_trace=False, record_stats=False)
            after = time.clock()
            times[k] += after - before
            print '\t', after-before, times[k]

    #return 0.05, np.ones(num_ints) 
    return hs, times

### SIMULATION PARAMETERS
#poly_name = 'tetrahedron'
#poly_name = 'octahedron'
poly_name = 'icosahedron'


verts, face_inds, cents = getattr(poly, poly_name)()
V, E, F, S, species, f_types, adj_list, dual = bga.get_poly(poly_name)
ints, ids, paths, shell_int, shell_paths, edges, shell_edge = bga.get_bg_ss(poly_name)

num_processes = 4

err_tol = 10**-9
N = 10**7
archive_rate = 10**6
output_rate = 10**5
h = 0.05

sample_time = 1.0
num_ints = len(ints)
int_list = range(1, num_ints)

run_str = "icosA"

if __name__ == "__main__":
    print_sysinfo()
    print get_time_str(), '\n'

    # Get sample rates, find trivial intermediates
    #h, sample_rates = optimize_sampling(poly_name, num_ints, err_tol)
    #hs, times = optimize_sampling(poly_name, num_ints, err_tol)
    sample_rates, skip_sampling = get_sample_rates(poly_name, int_list, h,err_tol, 10**3)
    
    #for int_num in range(1,len(ints)):
    #    sample_rates[int_num] = get_sample_rate(int_num, h, sample_time)

    # Split up ints
    int_groups = group_ints(num_processes, np.arange(1,len(ints)), weights=sample_rates)
        
    processes = [mp.Process(target=run_ints, 
                            args=(poly_name,
                                  int_groups[k], 
                                  N, 
                                  archive_rate, 
                                  output_rate, 
                                  run_str, 
                                  k,
                                  h,
                                  err_tol,
                                  skip_sampling)) for k in range(num_processes)]

    for p in processes:
        p.start()
    for p in processes:
        p.join()
    

