import numpy as np
from scipy.spatial.distance import pdist, squareform

import matplotlib.pyplot as plt
import scipy.integrate as integrate
import matplotlib.animation as animation

class Unit:
    def __init__(self, name, x, y, size):
        self.x = x
        self.y = y
        self.name = name
        self.size = size


#------------------------------------------------------------
# set up initial state


#------------------------------------------------------------
# set up figure and animation
fig = plt.figure()
fig.subplots_adjust(left=0, right=1, bottom=0, top=1)
ax = fig.add_subplot(111, aspect='equal', autoscale_on=False,
                     xlim=(0, 200), ylim=(0, 200))

# particles holds the locations of the particles
units1, = ax.plot([], [], 'bo', ms=6)
units2, = ax.plot([], [], 'ro', ms=6)

# rect is the box edge
bounds = [0,200,0,200]
rect = plt.Rectangle(bounds[::2],
                     bounds[1] - bounds[0],
                     bounds[3] - bounds[2],
                     ec='none', lw=2, fc='none')
ax.add_patch(rect)

def init():
    """initialize animation"""
    global rect
    units1.set_data([], [])
    units2.set_data([], [])
    rect.set_edgecolor('none')
    return units1, units2, rect

def animate(i):
    """perform animation step"""
    global rect, ax, fig

    ms = 6
    # update pieces of the animation
    rect.set_edgecolor('k')
    units1.set_data([10,50], [30,50])
    units2.set_data([110,150], [130,150])
    units1.set_markersize(ms)
    units2.set_markersize(ms)
    return units1, units2, rect

ani = animation.FuncAnimation(fig, animate, frames=600,
                              interval=10, blit=True, init_func=init)


# save the animation as an mp4.  This requires ffmpeg or mencoder to be
# installed.  The extra_args ensure that the x264 codec is used, so that
# the video can be embedded in html5.  You may need to adjust this for
# your system: for more information, see
# http://matplotlib.sourceforge.net/api/animation_api.html
ani.save('simulation.mp4', fps=30, extra_args=['-vcodec', 'libx264'])

#plt.show()
