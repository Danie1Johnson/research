import numpy as np
import matplotlib.pyplot as plt
import numpy.linalg
import cPickle
import os
import sys

#import manifold_reflected_brownian_motion as mrbm
import manifolds as mfs
import boundaries as bds
import statistics as sts
import bga_4_0 as bga
import multi_hist as mh

#mrbm = reload(mrbm)
mfs = reload(mfs)
bds = reload(bds)
sts = reload(sts)
mh = reload(mh)
#bga = reload(bga)


class MRBM:
    """
    Reflected Brownian Motion

    Simulate a RBM via a Random walk on a manifold.
    """
    
    def __init__(self, 
                 x0, 
                 h, 
                 manifold_name=None, 
                 unary_boundary_name=None, 
                 binary_boundary_name=None, 
                 stat_name=None,
                 manifold_kwargs={}, 
                 unary_boundary_kwargs={}, 
                 binary_boundary_kwargs={}, 
                 stat_kwargs={}, 
                 record_hist=False,
                 hist_min=None,
                 hist_max=None, 
                 hist_bins=None,
                 err_tol=10**-15, 
                 Sig=None):
        """
        """
        self.err_tol = err_tol
        

        # Manifold and Boundary functions
        self.c, self.C, self.manifold_reframe, self.manifold_mod_directions = mfs.get_manifold(manifold_name, 
                                                                                               kwargs=manifold_kwargs)

        self.unary_boundary = bds.get_boundary(unary_boundary_name, kwargs=unary_boundary_kwargs)      
        self.binary_boundary = bds.get_boundary(binary_boundary_name, kwargs=binary_boundary_kwargs, binary=True) 
        self.manifold_name = manifold_name
        self.unary_boundary_name = unary_boundary_name
        self.binary_boundary_name = binary_boundary_name
        self.manifold_kwargs = manifold_kwargs
        self.unary_boundary_kwargs = unary_boundary_kwargs
        self.binary_boundary_kwargs = binary_boundary_kwargs

        # Verify initial point is valid
        x0 = self.recenter(x0)
        if np.linalg.norm(self.c(x0)) > err_tol:
            raise Exception("ERROR: initial point x0 does not satisfy constraints ("
                            + str(numpy.linalg.norm(self.c(x0)))
                            + ") to tolerance ("
                            + str(self.err_tol)
                            + ").")
        if self.unary_boundary(x0) == False:
            raise Exception("ERROR: initial point x0 outside of unary_boundary")
        
        # Statistic functions & variables
        self.stat_name = stat_name
        if stat_name != None:
            #self.log_stats = log_stats
            self.stat, self.stat_strings = sts.get_stat(stat_name, kwargs=stat_kwargs)
            self.stat_sum = np.zeros_like(self.stat(x0))
            self.stat_log = np.array([self.stat(x0)])
            self.num_stats = len(self.stat(x0))
            self.record_hist = record_hist
            if record_hist == True:
                self.hist_min = hist_min
                self.hist_max = hist_max
                self.hist_bins = hist_bins
                self.hist = mh.MultiHist(self.hist_min, self.hist_max, self.num_stats, self.hist_bins)
    
        # Variables
        self.n = len(x0)
        #self.m = self.n - self.C(x0).shape[0]
        self.m = self.n - self.C(x0).shape[0] ############################################
        #print 'n:', self.n, '\tm:', self.m, self.C(x0).shape
        self.x0 = np.copy(x0)
        self.x = np.copy(x0)
        self.xs = np.array([np.copy(x0)])
        self.h = h
        self.T = np.array([0.0])
        self.samples = 1
        self.d = len(x0)
        if Sig == None:
            self.Sig = np.eye(self.m)
        else:
            self.Sig = Sig
        self.Sig_inv = numpy.linalg.inv(self.Sig)
        self.B = numpy.random.uniform(size=(self.n,self.m))

        #### Handle BG vs not, better ### need to access bg info regularly? or just once?
        #if self.unary_boundary_name == 'building_game':
        #    self.q0, self.links, self.lengths, self.faces = bga.load_bg_int(manifold_kwargs['poly_name'], 
        #manifold_kwargs['int_num'])

    def sample(self, 
               N=1, 
               record_trace=True, 
               record_stats=True,
               progress_bar=True, 
               bar_width=20):
        """
        Take a time grid and boundary axis lengths and draw samples according to rejection scheme.
        """
        T_run = numpy.arange(0.0,(N + 0.5)*self.h**2, self.h**2)
        
        if record_trace == True:
            xs_run = np.zeros((len(T_run), self.n))
            xs_run[0,:] = self.x
        if record_stats == True:
            stat_log_run = np.zeros((len(T_run), self.num_stats))
            stat_log_run[0,:] = self.stat(self.x)

        if progress_bar == True:
            #sys.stdout.write("[%s]" % (" " * bar_width))
            sys.stdout.write("_" * bar_width)
            sys.stdout.write('\n')
            sys.stdout.flush()
            #sys.stdout.write("\b" * (bar_width+1)) 
        else:
            bar_width = 1

        dN = N/(bar_width - 1)
        rN = N%(bar_width - 1)
        for b in xrange(bar_width):
            k_min = 1+dN*b
            if b == bar_width - 1:
                k_max = k_min + rN
            else:
                k_max = k_min + dN
            #for kt, t in enumerate(T_run[k_min:k_max]):
            for kt in xrange(k_min,k_max):
                self.x = self.new_rejection_sample()
                self.samples += 1        
                if record_trace == True:
                    xs_run[kt,:] = self.x
                if self.stat != None:
                    curr_stat = self.stat(self.x)
                    self.stat_sum += curr_stat
                    if record_stats == True:
                        stat_log_run[kt,:] = self.stat(self.x)
                    if self.record_hist == True:
                        self.hist.add_stats(curr_stat)
            if progress_bar == True:
                sys.stdout.write("-")
                sys.stdout.flush()
        if progress_bar == True:
            sys.stdout.write("\n")
        
        
        if record_trace == True:
            self.T = np.hstack((self.T, self.T[-1] + T_run[1:]))
            self.xs = np.vstack((self.xs, xs_run[1:,:]))

        if record_stats == True:
            self.stat_log = np.vstack((self.stat_log, stat_log_run[1:,:]))
        
        self.stat_mean = self.stat_sum/float(self.samples)

        return self.x

    def new_rejection_sample(self):
        """
        Draw sample according to rejection scheme.
        """
        x = np.copy(self.x)
        # Find Bases
        if self.manifold_mod_directions != None:
            D = self.manifold_mod_directions(x)
            rankD = D.shape[0]
            A = np.hstack((self.C(x).T, D.T, self.B[:,rankD:]))
        else:
            A = np.hstack((self.C(x).T, self.B))
            rankD = 0
        Q, R = numpy.linalg.qr(A)
        # Check for full rank
        if Q.shape[1] != self.n:
            raise Exception("ERROR?: A not of full rank", n, Q.shape[1])  
             
        Q1 = np.copy(Q[:,:self.n-self.m])
        Q1b = np.copy(Q[:,:self.n-self.m+rankD])
        Q2 = np.copy(Q[:,self.n-self.m:])
        Q2b = np.copy(Q[:,self.n-self.m+rankD:])

        x_prop = None
        while x_prop == None:
            gamma_sol = None
            while gamma_sol == None:
                # Make step in tangent space.
                alpha = numpy.random.multivariate_normal(np.zeros(self.m - rankD), np.eye(self.m - rankD))
                v = np.dot(Q2b,alpha)
                v /= numpy.linalg.norm(v)
                y = x + self.d**0.5*self.h*v 
                #print v
                # Project back to M
                #gamma = np.zeros(self.n-self.m + rankD)
                #F = lambda gam: self.c(y + np.dot(Q1b,gam))
                #J = lambda gam: np.dot(self.C(y + np.dot(Q1b,gam)),Q1b)
                #print gamma.shape, Q1b.shape, y.shape
                #print D.shape, rankD
                #print self.m, self.n
                #print self.c(y).shape, self.C(y).shape
                #print F(gamma).shape, J(gamma).shape
                gamma = np.zeros(self.n-self.m)
                F = lambda gam: self.c(y + np.dot(Q1,gam))
                J = lambda gam: np.dot(self.C(y + np.dot(Q1,gam)),Q1)
                gamma_sol = Newton(gamma, F, J, self.err_tol)
            #x_prop = y + np.dot(Q1b,gamma_sol)
            x_prop = y + np.dot(Q1,gamma_sol)

            # Check if x_prop within boundaries.
            if self.unary_boundary(x_prop) == False or self.binary_boundary(self.x, x_prop) == False:
                x_prop = None

        # If rotation controlled in the manifold, rotate back to fram of reference. 
        if self.manifold_reframe != None:####################################
            x_prop = self.manifold_reframe(x_prop, self.x)

        return x_prop

    def parm_str(self, N=None, M=None, t=None):
        pstr = self.manifold_name + "_" 
        pstr += "n_" + str(self.n) + "_"
        pstr += "m_" + str(self.m) + "_"
        for a in self.manifold_kwargs.items():
            pstr += a[0] + "_" + parm_to_str(a[1]) + "_"
        pstr += self.unary_boundary_name + "_"  
        for a in self.unary_boundary_kwargs.items():
            pstr += a[0] + "_" + parm_to_str(a[1]) + "_"
        pstr += "h_" + parm_to_str(self.h) + "_"

        if N != None:
            pstr += "N_" + str(N)
        if M != None:
            pstr += "M_" + str(M)
        if t != None:
            pstr += "t_" + parm_to_str(t)
        
        return pstr[:-1]

    def dump_trace(self):
        filename = self.pkl_filename()
        if os.path.isfile(filename):
            print "File", filename, "already exists. Nothing written."
        else:
            cPickle.dump(self.xs, open(filename, 'wb'))

    def load_trace(self):
        filename = self.pkl_filename()
        try:
            self.xs = cPickle.load(open(filename, 'rb'))
            return True
        except IOError:
            return False

    def pkl_filename(self):
        return self.parm_str(**self.run_args) + ".pkl"

    def recenter(self, x):
        if 'fixed_com' in self.manifold_kwargs:
            if self.manifold_kwargs['fixed_com'] == True:
                if 'masses' in self.manifold_kwargs:
                    masses = self.manifold_kwargs['masses']
                else:
                    masses = np.ones(len(x)/3)
                for j in range(3):
                    dim_com = np.dot(x[j::3], masses)/sum(masses)
                    x[j::3] -= dim_com
        return x


    def write_histograms(self, filename):
        """
        Write recorded histograms to file 
        """
        with open(filename, 'wb') as f:
            for k in range(self.hist.num_stats):
                if self.stat_strings != None:
                    f.write(self.stat_strings[k] + '\n')
                f.write(self.hist.hist_str(k) + '\n')
        


def project_to_hyperplane(norm_amb, z, Q, x):
    """
    Project the normal vector in the ambient space z_norm_amb at z_min 
    onto the plane tangent to the manifold at x with orthonormal basis Q.
    """
    n, m = Q.shape
    w = np.zeros(n)
    for k in range(m):
        w += np.dot(Q[:,k],z + norm_amb - x)*Q[:,k] 
    n_hat = x + w - z
    n_hat /= numpy.linalg.norm(n_hat)

    return n_hat

def Newton(x0, F, J, err_tol, max_itr=100):
    """
    Perform multidimensional Newton-Raphson iteration.
    """
    x = np.copy(x0)

    for k in xrange(max_itr):
        #print x.shape
        Fx = F(x)
        #print J(x), Fx
        if numpy.linalg.norm(Fx) < err_tol:
            return x
        dx = numpy.linalg.solve(J(x),-Fx)
        x += dx
        #yprint k
        
    raise Exception("ERROR: Maximum iteration exceeded in 'Newton'.")
        
def parm_to_str(x, dec_len=4):
    if isinstance(x, np.ndarray):
        pstr = ""
        for y in x:
            pstr += parm_to_str(y, dec_len) + "_"
        return pstr[:-1]
    else:
        if x < 0:
            print "ERROR: negative parameter."
            return None
        elif x >= 10.0:
            print "ERROR: parameters > 10 not valid."
            return None

        x_str = str(x)    
        return x_str[0]+"_"+x_str[2:min(2+dec_len,len(x_str))]


