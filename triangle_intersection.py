import numpy as np

def check_intersection(t1, t2):
    """
    Take two length 3x3 numpy arrays with rows representing the 3 vertices. 
    Return True if they intersect and False otherwise. 
    Use Algorithm in 'A Fast Triangle-Triangle Intersetion Test' by Tomas Moller.
    """
    
    # 1. Compute t2 plane.
    N2 = np.cross(t2[1,:] - t2[0,:],t2[2,:] - t2[0,:])

    # 2. Reject if all points of t1 on same side.
    d10 = np.dot(N2,t1[0,:] - t2[0,:])
    d11 = np.dot(N2,t1[1,:] - t2[0,:])
    d12 = np.dot(N2,t1[2,:] - t2[0,:])
    if d10*d11 >= 0 and d10*d12 >= 0:
        return False

    # 3. Compute t1 plane.
    N1 = np.cross(t1[1,:] - t1[0,:],t1[2,:] - t1[0,:])

    # 4. Reject if all points of t2 on same side.
    d20 = np.dot(N1,t2[0,:] - t1[0,:])
    d21 = np.dot(N1,t2[1,:] - t1[0,:])
    d22 = np.dot(N1,t2[2,:] - t1[0,:])
    if d20*d21 >= 0 and d20*d22 >= 0:
        return False

    # 5. Compute intersection line and project onto largest axis.
    D = np.cross(N1,N2)
    D_abs = abs(D)

    if D_abs[0] > D_abs[1] and D_abs[0] > D_abs[2]:
        p10 = t1[0,0]
        p11 = t1[1,0]
        p12 = t1[2,0]
    elif D_abs[1] > D_abs[2]:
        p10 = t1[0,1]
        p11 = t1[1,1]
        p12 = t1[2,1]
    else:
        p10 = t1[0,2]
        p11 = t1[1,2]
        p12 = t1[2,2]

    # 6. Compute intervals for each triangle.
    x1 = p10 + ( 
    # 7. Intersect the intervals 
    
