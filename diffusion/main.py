
# coding: utf-8

# In[34]:

# import numpy as np
# import math as m
# import matplotlib.pyplot as plt
# from scipy.integrate import odeint
# from functools import reduce
# import csv
# # %matplotlib inline
# from scipy.interpolate import *


import sys
import os

sys.path.append(os.path.abspath('..'))


from utils import *

get_ipython().magic('matplotlib inline')
from diffussion_template import *
from animation_helper import *
from matplotlib import rc
# plt.rcParams['animation.ffmpeg_path']
os.path.abspath('../../ffmpeg')
# sys.path

plt.rcParams['animation.ffmpeg_path'] = os.path.abspath('../../ffmpeg/bin')




# In[35]:

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

# Variable / Parameter names
#
# n_cells #-- size of grid
# x_coords = np.arange(0,n_cells) #-- x_coordinates of each cell
#
# grid #-- array of values representing concentration in each cell
# kD #-- diffusion rate
#
def d_diffusion(grid,kD):

    n_cells=len(grid)
    d_diff = np.zeros(n_cells)

    # d_diff[i] is the rate of change of
    # conc. in cell grid[i]

    # the rate of flow out of each cell
    # into each neighbour
    # is equal to kD*[conc in cell]

    # we need to use a loop over each cell
    # in the 1D grid and calculate the net rate of
    # change due to diffusion

    # hint
    # 1. write a loop over each cell
    # 2. identify the neighbouring cells
    # 3. sum the rates of flows in and out of this cell
    # 4. store this value in the appropriate element of
    #     array d_diff
#     d_grad=n_cells*[0];
    c0=grid[-1];
    c1=grid[0];
    for i in range(n_cells):
        c2=grid[(i+1)%n_cells];
        d_diff[i]= kD*( -2.*c1+c0+c2);
        c0=c1;
        c1=c2;
    

#     **DELETE THIS AND INSERT CODE FOR DIFFUSION**

    return d_diff

def sdot(s,t,params):

    (kD,)=params
    # kD is the rate of diffusion

    # s is an array of values
    # for each cell in the grid

    # find size of grid
    n_cells=len(s)

    # create rate array for results
    ds=np.zeros(n_cells)

    # process any reactions occuring in each cell
    for i in range(n_cells):
        pass # the pass command allows us to have an empty loop
             # meaning 'do nothing' this can
             # act as a placeholder for when we want to
             # add reactions e.g. ds[i]= -k*s[i] for degradation

    # calc net rate diffusion in each cell
    # using external function that acts on grid s
    ds=ds+d_diffusion(s,kD)

    return ds






# In[36]:

n_cells=20
x_coords=np.arange(0,n_cells)

# create an array for initial cell conditions
s0=np.zeros(n_cells)
# initialises grid with each cell conc set to 0.0

# **DELETE THIS BUT SET CORRECT INITIAL CONDITIONS HERE**
s0=np.ones(n_cells);
s0[4]=3.0;
# s0=np.ones(n_cells)

# set all initial concentrations to 1.0
# except in cell 5 which has concentration 3.0


# set diffusion rate
# same for all cells
kD=0.2

# create array of time observations
t_min=0
t_max=800.
t_interval=10.
t_obs=np.arange(t_min,t_max,t_interval)

# set initial conditions and parameters
s=s0
params=(kD,)

# run simulation
s_obs=odeint(sdot,s,t_obs,args=(params,),mxstep=5000000)
# mxstep=5000000 allows more calculations to be taken per timestep
# otherwise Python may complain it is taking too long to solve

# convert to array so we can manipulate it more easily
s_obs=np.array(s_obs)

# plot results
plt.close("all")
fig1 = plt.figure()
ax=plt.subplot(1,2,1)

# each row of s_obs is a time observation
# containing a column for each cell
n_obs=len(t_obs)
for i in range(n_obs):
    ax.plot(x_coords, s_obs[i,:])

ax.set_ylim([0.8,3.3])
fig1.show()

fig2=plt.figure();
ax2=fig1.add_subplot(1,2,2)
line,=ax2.plot(x_coords,s_obs[0,:]);
ax2.set_ylim([0.8,3.2])
anim=plot_time_series(fig2,line,s_obs,frame_interval=5);
fig2.show()


# In[43]:



from tempfile import NamedTemporaryFile

VIDEO_TAG = """<video controls>
 <source src="data:video/x-m4v;base64,{0}" type="video/mp4">
 Your browser does not support the video tag.
</video>"""

def anim_to_html(anim):
    if not hasattr(anim, '_encoded_video'):
        with NamedTemporaryFile(suffix='.mp4') as f:
            anim.save(f.name, fps=20,  writer = FFwriter, extra_args=['-vcodec', 'libx264'])
            video = open(f.name, "rb").read()
        anim._encoded_video = video.encode("base64")
    
    return VIDEO_TAG.format(anim._encoded_video)


from IPython.display import HTML

def display_animation(anim):
    plt.close(anim._fig)
    return HTML(anim_to_html(anim))

#FFwriter = animation.FFMpegWriter()
#display_animation(anim)
#exce(tst)
# from IPython.display import HTML


# 'mencoder'
# anim._repr_html_() is None
# rc('animation', html='html5')

# HTML(anim.to_html5_video())
# anim


# FFwriter.__dict__

