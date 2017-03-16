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
    d_grad=n_cells*[0];
    c0=grid[-1];
    c1=grid[0];
    for i in range(grid):
        c2=grid[(i+1)%n_cells];
        d_grad[i]= kD*( -2.*c1+c0+c20);
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



