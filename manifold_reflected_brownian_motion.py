import numpy as np
import matplotlib.pyplot as plt
import numpy.linalg
import cPickle
import os

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
                 manifold_name, 
                 boundary_name, 
                 x0, 
                 h, 
                 stat_name=None,
                 record_hist=False,
                 hist_min=None, 
                 hist_max=None, 
                 hist_bins=None,
                 manifold_kwargs={}, 
                 boundary_kwargs={}, 
                 stat_kwargs={}, 
                 err_tol=10**-15, 
                 Sig=None):
        """
        """
        # Manifold and Boundary functions
        self.c, self.C = mfs.get_manifold(manifold_name, kwargs=manifold_kwargs)
        self.boundary = bds.get_boundary(boundary_name, kwargs=boundary_kwargs)      
        self.manifold_name = manifold_name
        self.boundary_name = boundary_name
        self.manifold_kwargs = manifold_kwargs
        self.boundary_kwargs = boundary_kwargs
        
        # Statistic functions & variables
        self.stat_name = stat_name
        if stat_name != None:
            #self.log_stats = log_stats
            self.stat = sts.get_stat(stat_name, kwargs=stat_kwargs)
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
        self.m = self.n - self.C(x0).shape[0]
        self.x0 = np.copy(x0)
        self.x = np.copy(x0)
        self.xs = np.array([np.copy(x0)])
        self.h = h
        self.T = np.array([0.0])
        self.err_tol = err_tol
        self.samples = 1
        self.d = len(x0)
        if Sig == None:
            self.Sig = np.eye(self.m)
        else:
            self.Sig = Sig
        self.Sig_inv = numpy.linalg.inv(self.Sig)
        self.B = numpy.random.uniform(size=(self.n,self.m))

        ### Handle BG vs not, better 
        if self.boundary_name != 'none':
            self.dihedral_inds = []
            self.q0, self.links, self.lengths, self.faces = bga.load_bg_int(manifold_name, manifold_kwargs['int_num'])
            F = len(self.faces)
            for j in range(F):
                for k in range(j+1, F):
                    if bds.adjacent_faces(j, k, self.faces) == True:
                        self.dihedral_inds.append(self.order_verts(self.faces[j], self.faces[k]))
            self.dihedrals = self.get_dihedrals(x0)

        if self.boundary(x0) == False:
            print "ERROR: x0 not in domain." 


    def sample(self, N=1, record_trace=True, record_stats=True):
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

        for kt, t in enumerate(T_run[1:]):
            self.x = self.new_rejection_sample()
            if record_trace == True:
                xs_run[kt+1,:] = self.x
            if self.stat != None:
                curr_stat = self.stat(self.x)
                self.stat_sum += curr_stat
                if record_stats == True:
                    stat_log_run[kt+1,:] = self.stat(self.x)
                if self.record_hist == True:
                    self.hist.add_stats(curr_stat)
            if self.boundary_name != 'none':
                self.dihedrals = self.get_dihedrals(self.x)

        self.samples += N        

        if record_trace == True:
            self.T = np.hstack((self.T, self.T[-1] + T_run[1:]))
            self.xs = np.vstack((self.xs, xs_run[1:,:]))

        if record_stats == True:
            self.stat_log = np.vstack((self.stat_log, stat_log_run[1:,:]))
        
        return self.x

    def new_rejection_sample(self):
        """
        Draw sample according to rejection scheme.
        """
        x = np.copy(self.x)
        # Find Bases
        A = np.hstack((self.C(x).T, self.B))
        Q, R = numpy.linalg.qr(A)
        # Check for full rank
        if Q.shape[1] != self.n:
            raise Exception("ERROR?: A not of full rank", n, Q.shape[1])  
             
        Q1 = np.copy(Q[:,:self.n-self.m])
        Q2 = np.copy(Q[:,self.n-self.m:])

        x_prop = None
        while x_prop == None:
            gamma_sol = None
            while gamma_sol == None:
                # Make step in tangent space.
                alpha = numpy.random.multivariate_normal(np.zeros(self.m),self.Sig)
                v = np.dot(Q2,alpha)
                v /= numpy.linalg.norm(v)
                y = x + self.d**0.5*self.h*v 

                # Project back to M
                gamma = np.zeros(self.n-self.m)
                F = lambda gam: self.c(y + np.dot(Q1,gam))
                J = lambda gam: np.dot(self.C(y + np.dot(Q1,gam)),Q1)
                gamma_sol = Newton(gamma, F, J, self.err_tol)
            x_prop = y + np.dot(Q1,gamma_sol)

            # Check if x_prop within boundary.
            if self.boundary(x_prop) == False or (self.boundary_name != 'none' and 
                                                  self.dihedral_switch(x_prop) == True):
                x_prop = None

        return x_prop

    def order_verts(self, v1, v2):
        """
        get order of vert inds for dihedral calculation.
        """
        common_verts = set(v1).intersection(set(v2))
        ordered_verts = [list(set(v1).difference(common_verts))[0],
                         list(common_verts)[0],
                         list(common_verts)[1],
                         list(set(v2).difference(common_verts))[0]]
        return ordered_verts

    def dihedral_switch(self, y):
        """
        Compare the dihedrals of proposal point y with self.x and check for sign change. 
        """
        new_dihedrals = self.get_dihedrals(y)

        #if min((self.dihedrals - np.pi)*(new_dihedrals - np.pi)) < -1.0:
        
        if max(abs(new_dihedrals - self.dihedrals)) > 2.0:
            return True
            #return False ### NO NO NO 
        else:
            return False

    def get_dihedrals(self, y):
        """
        Compure the dihedral angles at y.
        """
        dihedrals = []
        for di in self.dihedral_inds:
            dihedrals.append(sts.signed_dihedral_angle(y, di[0], di[1], di[2], di[3]))
        return np.array(dihedrals)

    def parm_str(self, N=None, M=None, t=None):
        pstr = self.manifold_name + "_" 
        pstr += "n_" + str(self.n) + "_"
        pstr += "m_" + str(self.m) + "_"
        for a in self.manifold_kwargs.items():
            pstr += a[0] + "_" + parm_to_str(a[1]) + "_"
        pstr += self.boundary_name + "_"  
        for a in self.boundary_kwargs.items():
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


