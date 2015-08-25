import numpy as np
from scipy.spatial.distance import pdist, squareform

import matplotlib.pyplot as plt
import scipy.integrate as integrate
import matplotlib.animation as animation
import math


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
        if str and str != "-":
            tmp = str.split(",")
            data_x.append(float(tmp[0]))
            data_y.append(float(tmp[1]))
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
                     xlim=(0, 200), ylim=(0, 200))



len = (file_len('./paths1.txt')-2)#/10



file1 = open('./paths1.txt', 'r')
file2 = open('./paths2.txt', 'r')

#units1 = getUnitStats(file1)
#units2 = getUnitStats(file2)

sizes1 = getUnitStats(file1)
sizes2 = getUnitStats(file2)

#for unit in units1:
#    sizes1.append(math.sqrt(unit.size))
#for unit in units2:
#    sizes2.append(math.sqrt(unit.size))

positions1 = []
positions2 = []



for s in sizes1:
    tmp, = ax.plot([],[],'ro', ms=s*3)
    positions1.append(tmp)
for s in sizes2:
    tmp, = ax.plot([],[],'bo', ms=s*3)
    positions2.append(tmp)
    


# rect is the box edge
bounds = [0,200*5,0,200*5]
rect = plt.Rectangle(bounds[::2],
                     bounds[1] - bounds[0],
                     bounds[3] - bounds[2],
                     ec='none', lw=2, fc='none')
ax.add_patch(rect)

def init():
    """initialize animation"""
    global rect
    for pos in positions1:
        pos.set_data([],[])
    for pos in positions2:
        pos.set_data([],[])
    
    rect.set_edgecolor('none')
    return positions1, positions2, rect

def animate(i):
    """perform animation step"""
    global rect, ax, fig, file1, file2, units1, units2

    # update pieces of the animation
    rect.set_edgecolor('k')
    
        
#    for i in xrange(0,9):
#        file1.readline()
#        file2.readline()

    a,b = getData(file1)
    for pos, data in zip(positions1, zip(a,b)):
        pos.set_data(data[0],data[1])
    a,b = getData(file2)
    for pos, data in zip(positions2, zip(a,b)):
        pos.set_data(data[0],data[1])
        
    return positions1, positions2, rect

ani = animation.FuncAnimation(fig, animate, frames=len,
                              interval=10, blit=True, init_func=init)


# save the animation as an mp4.  This requires ffmpeg or mencoder to be
# installed.  The extra_args ensure that the x264 codec is used, so that
# the video can be embedded in html5.  You may need to adjust this for
# your system: for more information, see
# http://matplotlib.sourceforge.net/api/animation_api.html
ani.save('simulation.mp4', fps=30, extra_args=['-vcodec', 'libx264'])

file1.close()
file2.close()

#plt.show()
