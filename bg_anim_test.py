import numpy as np
import matplotlib.pyplot as plt
from time import time
import matplotlib.animation as animation

import bga_4_0 as bga
import manifold_reflected_brownian_motion as mrbm

bga = reload(bga)
mrbm = reload(mrbm)

poly_name = 'octahedron'
int_num = 5
boundary_name = 'none'
manifold_name = poly_name + "__" + str(int_num)

scheme = 'rej'
h = 0.01
N = 10

n, dim, x0, masses, links, lengths, faces = bga.load_bg_int(poly_name, int_num)
z = mrbm.MRBM(manifold_name, boundary_name, x0, h, scheme, run_args={'N': N})
           
z.sample(N=1000)


#------------------------------------------------------------
# set up figure and animation
fig = plt.figure()
ax = fig.add_subplot(111, aspect='equal', autoscale_on=False,
                     xlim=(-2.0, 2.0), ylim=(-2.0, 2.0))
ax.grid()

line, = ax.plot([], [], 'o-', lw=2)
line2, = ax.plot([], [], 'o-', lw=2)
#line3, = ax.plot([], [], 'o-', lw=2)
#line4, = ax.plot([], [], 'o-', lw=2)
#time_text = ax.text(0.02, 0.95, '', transform=ax.transAxes)
#energy_text = ax.text(0.02, 0.90, '', transform=ax.transAxes)
#beta_text = ax.text(0.02, 0.85, '', transform=ax.transAxes)


def face_position(bg_int, face_num, faces, dim=3):
    """Return the current and last positions in the desired dimensions for viewing."""

    x_list = []
    y_list = []          
    ## Right now, only views in x and y dimensions
    for v in faces[face_num]:
        x_list.append(bg_int.x[dim*v])
        y_list.append(bg_int.x[dim*v + 1])

    x_list.append(bg_int.x[dim*faces[face_num][0]])
    y_list.append(bg_int.x[dim*faces[face_num][0] + 1])

    return (np.array(x_list), np.array(y_list))


def init():
    """initialize animation"""
    line.set_data([], [])
    line2.set_data([], [])
    #line3.set_data([], [])
    #line4.set_data([], [])
    #time_text.set_text('')
    #energy_text.set_text('')
    #beta_text.set_text('')
    return line, line2, #line3, time_text, energy_text, beta_text


def animate_bg(i, bg_int, faces):
    """perform animation step"""
    bg_int.sample()
    ## Fix for multi face animations.
    line.set_data(face_position(bg_int, 0, faces))
    line2.set_data(face_position(bg_int, 1, faces))
    #time_text.set_text('time = %.1f' % tle.time_elapsed)
    #energy_text.set_text('residual = %.6f' % tle.res())
    #beta_text.set_text('estimate = %.4f' % tle.angle_occupation(0, ang_val, ang_tol))
    return line, line2, #line3, line4, time_text, energy_text, beta_text


def bg_animation(fig, bg_int, faces, save_animation=False):
    """
    """
    animate = lambda x: animate_bg(x, bg_int, faces)
    L = 1.0

    face_lines = ()

    # choose the interval based on dt and the time to animate one step
    t0 = time()
    animate(0)
    t1 = time()
    interval = 1000 * (L * bg_int.h) - (t1 - t0)


    ani = animation.FuncAnimation(fig, animate, frames=300,
                                  interval=interval, blit=True, init_func=init)

    if save_animation == True:
        # save the animation as an mp4.  This requires ffmpeg or mencoder to be
        # installed.  The extra_args ensure that the x264 codec is used, so that
        # the video can be embedded in html5.  You may need to adjust this for
        # your system: for more information, see
        # http://matplotlib.sourceforge.net/api/animation_api.html
        ani.save('triangular_linkage_diffusion.mp4', fps=3, extra_args=['-vcodec', 'libx264'])

    plt.show()




bg_animation(fig, z, faces)
