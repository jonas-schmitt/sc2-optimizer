import numpy as np
from scipy.spatial.distance import pdist, squareform

import matplotlib.pyplot as plt
import scipy.integrate as integrate
import matplotlib.animation as animation
import math
import sys


class Unit:
    def __init__(self, name, size):
        self.name = name
        self.size = size

def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1

def getData(f):
    line = f.readline()
    line = line.strip('\n')
    str_list = line.split('\t\t\t')
    data_x = []
    data_y = []
    for str in str_list:
        if str:
            if str != "-":
                tmp = str.split(",")
                data_x.append(float(tmp[0]))
                data_y.append(float(tmp[1]))
            else:
                data_x.append(-100)
                data_y.append(-100)
    return data_x, data_y

def getUnitStats(f):
    line = f.readline()
    line = line.strip('\n')
    str_list = line.split('\t\t\t')
    units = []
    #for str in str_list:
    #    units.append(Unit(str, 1))
    line = f.readline()
    line = line.strip('\n')
    str_list = line.split('\t\t\t')
    for str in str_list:
        if str:
            units.append(float(str))
    return units




#------------------------------------------------------------
# set up initial state


#------------------------------------------------------------
# set up figure and animation
fig = plt.figure()
fig.subplots_adjust(left=0, right=1, bottom=0, top=1)
ax = fig.add_subplot(111, aspect='equal', autoscale_on=False,
                     xlim=(0, 100), ylim=(0, 100))

paths = [str(sys.argv[1]), str(sys.argv[2])]

len = (file_len(paths[0])-2)#/10




file1 = open(paths[0], 'r')
file2 = open(paths[1], 'r')

#units1 = getUnitStats(file1)
#units2 = getUnitStats(file2)

sizes1 = getUnitStats(file1)
sizes2 = getUnitStats(file2)

#for unit in units1:
#    sizes1.append(math.sqrt(unit.size))
#for unit in units2:
#    sizes2.append(math.sqrt(unit.size))

circles1 = []
circles2 = []


# rect is the box edge
bounds = [0,100,0,100]
rect = plt.Rectangle(bounds[::2],
                     bounds[1] - bounds[0],
                     bounds[3] - bounds[2],
                     ec='none', lw=2, fc='none')

size = fig.get_size_inches()*fig.dpi
print size

for s in sizes1:
    circle = plt.Circle((0,0),s,color='r',fill=True, clip_on = False)
    circles1.append(circle)
    ax.add_artist(circle)
for s in sizes2:
    circle = plt.Circle((0,0),s,color='b',fill=True, clip_on = False)
    circles2.append(circle)
    ax.add_artist(circle)



ax.add_patch(rect)

def init():
    """initialize animation"""
    global rect
#    for circle in circles1:
#        circle.center = 0,0
#    for circle in circles2:
#        circle.center = 0,0

    rect.set_edgecolor('none')
    return circles1, circles2, rect

def animate(i):
    """perform animation step"""
    global rect, ax, fig, file1, file2, units1, units2

    # update pieces of the animation
    rect.set_edgecolor('k')


#    for i in xrange(0,9):
#        file1.readline()
#        file2.readline()

    a,b = getData(file1)
    for circle, data in zip(circles1, zip(a,b)):
#        if not data[0] or not data[1]:
#            pos.set_data([-100],[-100])
#        else:
	    circle.center = data[0],data[1]
    a,b = getData(file2)
    for circle, data in zip(circles2, zip(a,b)):
#        if not data[0] or not data[1]:
#            pos.set_data([-100],[-100])
#        else:
            circle.center = data[0],data[1]

    return circles1, circles2, rect

ani = animation.FuncAnimation(fig, animate, frames=len,
                              interval=100, blit=True, init_func=init)


# save the animation as an mp4.  This requires ffmpeg or mencoder to be
# installed.  The extra_args ensure that the x264 codec is used, so that
# the video can be embedded in html5.  You may need to adjust this for
# your system: for more information, see
# http://matplotlib.sourceforge.net/api/animation_api.html
ani.save('simulation.mp4', fps=30, extra_args=['-vcodec', 'libx264'])

file1.close()
file2.close()

#plt.show()
