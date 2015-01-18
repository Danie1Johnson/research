import numpy as np
import matplotlib.pyplot as plt
import numpy.linalg
import cPickle
import os

#import manifold_reflected_brownian_motion as mrbm
import manifolds as mfs
import boundaries as bds
import statistics as sts

#mrbm = reload(mrbm)
mfs = reload(mfs)
bds = reload(bds)
sts = reload(sts)
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
                 scheme, 
                 stats_name=None, 
                 run_kwargs={}, 
                 manifold_kwargs={}, 
                 boundary_kwargs={}, 
                 stat_kwargs={}, 
                 err_tol=10**-15, 
                 Sig=None):
        """
        """
        # Manifold and Boundary functions
        self.c, self.C = mfs.get_manifold(manifold_name, kwargs=manifold_kwargs)
        #self.boundary, self.boundary_normal = bds.get_boundary(boundary_name, kwargs=boundary_kwargs)
        self.boundary = bds.get_boundary(boundary_name, kwargs=boundary_kwargs)
        
        self.manifold_name = manifold_name
        self.boundary_name = boundary_name
        self.run_args = run_kwargs
        self.manifold_kwargs = manifold_kwargs
        self.boundary_kwargs = boundary_kwargs
        
        # Statistic functions & variables
        self.stat_name = stat_name
        if stat_name != None:
            self.log_stats = log_stats
            self.stats = sts.get_stats(stat_name, kwargs=stat_kwargs)
            self.stat_sums = np.zeros_like(stats(x0))
            self.stats_log = np.array([self.stats(x0)])
            self.num_stats = len(self.stats_log)
    
        # Variables
        self.n = len(x0)
        #print x0
        #print self.n, self.C(x0).shape
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
        self.scheme = scheme

        ## Functions
        #self.boundary = boundary
        #self.boundary_normal = boundary_normal
        #self.c = c
        #self.C = C
        
        # Verify that x0 is in the domain
        #self.N_bounds = len(self.boundary(x0))
        #if np.all(self.boundary(x0)) == False:
        if self.boundary(x0) == False:
            print "ERROR: x0 not in domain." 

    def boundary_check(self, x):
        """
        For each f(x) < 0.0 constraint in domain, return list of bounds not satisfied.
        """

        return np.where(self.boundary(x) > 0)[0]
    

    #def find_boundary_intersection(self, x, xp, b_num, max_itrs=100):
    #    """
    #    Use bisection to find the point z at which the segment 
    #    x <-> xp intersects boundary # b_num. 
    #
    #    i.e. find c \in [0,1] such that d*xp + (1-d)*x is 
    #    within the boundary for d < c and outside for d > c.
    #
    #    Since we only need to find z within err_tol of the boundary, 
    #    make sure z is just within the boundary.
    #    """
    #    c_mx = 1.0
    #    c_mn = 0.0
    #
    #    if (b_num in self.boundary_check(xp)) == False:
    #        print "ERROR: xp not outside boundary", xp
    #        return None
    #    if b_num in self.boundary_check(x):
    #        print "ERROR: x outside boundary", x
    #        return None
    #        
    #    for its in range(max_itrs):
    #        if numpy.linalg.norm((c_mx - c_mn)*xp + (-c_mx + c_mn)*x) < self.err_tol:
    #            return c_mn*xp + (1.0 - c_mn)*x, c_mn
    #        else:
    #            c_ave = 0.5*(c_mn + c_mx)
    #            c_ave_viol = self.boundary_check(c_ave*xp + (1.0 - c_ave)*x)
    #            if b_num in c_ave_viol:
    #                c_mx = c_ave
    #            else:
    #                c_mn = c_ave
    #    print "ERROR: max_itrs exceeded in find_boundary_intersection",  
    #    print numpy.linalg.norm((c_mx - c_mn)*xp + (-c_mx + c_mn)*x), "greater than", self.err_tol
    #    return None
                
    

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
            stat_log_run[0,:] = self.stats(self.x)

        for kt, t in enumerate(T_run[1:]):
            self.x = self.new_rejection_sample()
            if record_trace == True:
                xs_run[kt+1,:] = self.x
            if record_stats == True:
                stat_log_run[kt+1,:] = self.x
        self.samples += N        

        if record_trace == True:
            self.T = np.hstack((self.T, self.T[-1] + T_run[1:]))
            self.xs = np.vstack((self.xs, xs_run[1:,:]))

        if stats != None:
            self.stats_sum += self.stats(self.x)
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
            print "ERROR?: A not of full rank", n, Q.shape[1]  
            return
        Q1 = np.copy(Q[:,:self.n-self.m])
        Q2 = np.copy(Q[:,self.n-self.m:])

        # Make step in tangent space.
        #print self.m, self.Sig
        alpha = numpy.random.multivariate_normal(np.zeros(self.m),self.Sig)
        v = np.dot(Q2,alpha)
        v /= numpy.linalg.norm(v)
        y = x + self.d**0.5*self.h*v 

        gamma_sol = None
        while gamma_sol == None:
            # Check if xp within boundary.
            while self.boundary(y) == False:
                alpha = numpy.random.multivariate_normal(np.zeros(self.m),self.Sig)
                v = np.dot(Q2,alpha)
                v /= numpy.linalg.norm(v)
                y = x + self.d**0.5*self.h*v 

            # Project back to M
            gamma = np.zeros(self.n-self.m)
            F = lambda gam: self.c(y + np.dot(Q1,gam))
            J = lambda gam: np.dot(self.C(y + np.dot(Q1,gam)),Q1)
            gamma_sol = Newton(gamma, F, J, self.err_tol)
        x = y + np.dot(Q1,gamma_sol)

        return x


    #def new_reflection_sample(self):
    #    """
    #    Draw sample according to rejection scheme.
    #    """
    #    x = np.copy(self.x)
    #    L = (self.d)**0.5*self.h
    #    # Find Bases
    #    A = np.hstack((self.C(x).T, self.B))
    #    Q, R = numpy.linalg.qr(A)
    #    # Check for full rank
    #    if Q.shape[1] != self.n:
    #        print "ERROR?: A not of full rank", n, Q.shape[1]  
    #        return
    #    Q1 = np.copy(Q[:,:self.n-self.m])
    #    Q2 = np.copy(Q[:,self.n-self.m:])
    #
    #    # Make step in tangent space.
    #    alpha = numpy.random.multivariate_normal(np.zeros(self.m),self.Sig)
    #    v = np.dot(Q2,alpha)
    #    v /= numpy.linalg.norm(v)
    #    y = x + L*v 
    #
    #    # Check if y lies outside the boundary.
    #    x_viol = self.boundary_check(y)
    #    while len(x_viol) > 0:
    #        #print 'here'
    #        boundary_crossings = [self.find_boundary_intersection(x, y, b_num) for b_num in x_viol]
    #        # Find first boundary crossed (smallest c value)
    #        c_min = boundary_crossings[0][1]
    #        z_min = boundary_crossings[0][0]
    #        b_num_min = x_viol[0]
    #        for j, (z, c_z) in enumerate(boundary_crossings):
    #            if c_z < c_min:
    #                c_min = c_z
    #                z_min = z
    #                b_num_min = x_viol[j]
    #        if len(self.boundary_check(z_min)) > 0:
    #            print "ERROR: z_min outside.", z_min, c_min, b_num_min
    #        L -= numpy.linalg.norm(x - z_min) 
    #        
    #        # Get boundary normal at z_min.
    #        z_norm_amb = self.boundary_normal(z_min, b_num_min)
    #
    #        # Project to the tangent plane.
    #        z_norm = project_to_hyperplane(z_norm_amb, z_min, Q2, self.x)
    #        z_norm /= numpy.linalg.norm(z_norm)
    #
    #        #print z_norm, z_min, x, z_norm_amb
    #        # Find unit reflection vector
    #        v_reflect = z_min - x - 2.0*np.dot(z_min - x, z_norm)*z_norm
    #        #print v_reflect, Q2
    #        v_reflect /= numpy.linalg.norm(v_reflect)
    #        #print v_reflect
    #        x = z_min
    #        y = z_min + L*v_reflect
    #        # Check if y lies outside the boundary. 
    #        x_viol = self.boundary_check(y)
    #        
    #    # Project back to M
    #    gamma = np.zeros(self.n-self.m)
    #    F = lambda gam: self.c(y + np.dot(Q1,gam))
    #    J = lambda gam: np.dot(self.C(y + np.dot(Q1,gam)),Q1)
    #    gamma_sol = Newton(gamma, F, J, err_tol=self.err_tol)
    #    x = y + np.dot(Q1,gamma_sol)
    #
    #    if np.all(self.boundary(x) < 0.0) == False:
    #        print "ERROR: projection lies outside domain", x
    #        return None
    #
    #    return x

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
        pstr += self.scheme + "_"
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
    #print 'proj', norm_amb_z, Q
    n, m = Q.shape
    w = np.zeros(n)
    for k in range(m):
        w += np.dot(Q[:,k],z + norm_amb - x)*Q[:,k] 
    #print w
    n_hat = x + w - z
    #print 
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
        
    print "ERROR: Maximum iteration exceeded in 'Newton'."
    return
        
def rw_step(x, c, C, n, m, h, Sig, B, err_tol):
    """
    Take a random step of size h on Manifold defined by c and C with metric Sig.
    """
    # Find Bases
    #print C(x).shape, B.shape
    A = np.hstack((C(x).T, B))
    Q, R = numpy.linalg.qr(A)
    # Check for full rank
    if Q.shape[1] != n:
        print "ERROR?: A not of full rank", n, Q.shape[1]  
        return
    Q1 = np.copy(Q[:,:n-m])
    Q2 = np.copy(Q[:,n-m:])
    
    # Make step in tangent space.
    #if m > 1:
    alpha = numpy.random.multivariate_normal(np.zeros(m),Sig)
    #print alpha.shape, 222
    #if m == 1:

    #print Q2.shape, alpha.shape
    v = np.dot(Q2,alpha)
    v /= numpy.linalg.norm(v)
    y = x + h*v 

    # Project back to M
    gamma = np.zeros(n-m)
    
    #print Q1.shape, y.shape
    F = lambda gam: c(y + np.dot(Q1,gam))
    J = lambda gam: np.dot(C(y + np.dot(Q1,gam)),Q1)
    #print x.shape, 1
    gamma_sol = Newton(gamma, F, J, err_tol=err_tol)
    x = y + np.dot(Q1,gamma_sol)
    return x


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
