# -*- coding: utf-8 -*-
from matplotlib import pyplot as plt
from matplotlib.path import Path
import matplotlib.patches as patches
import numpy as np
import math

#####################################################
#
# INFORMATION
#
# The code below allows us to integrate
# dy/dt = t^3 - 1.5t^2 + 4
#
# It also display the results on a graph as the area
# under the integrated curve
#
######################################################


def dydt(t):
	return  t**3 - 3*t**2/2 + 4.

def y(t):
	return (t**4)/4 - (t**3)/2 + 4*t + 20

def area_of_strip_under_dydt(t_start,width):
	area=dydt(t_start)*width
	return area

def area_of_trapezium_under_dydt(t_start,width):
	t_end=t_start+width
	area=0.5 * ( dydt(t_start)+dydt(t_end) )*width
	return area

def sum_strips(t_start, t_end, n_strips):
	interval=t_end-t_start
	width=interval/n_strips
	total_area=0
	for i in range(n_strips):
		strip_start=t_start+ i *width
		total_area+=area_of_strip_under_dydt(strip_start,width)
	return total_area

def sum_trapezium_strips(t_start, t_end, n_strips):
	interval=t_end-t_start
	width=interval/n_strips
	total_area=0
	for i in range(n_strips):
		strip_start=t_start+ i *width
		total_area+=area_of_trapezium_under_dydt(strip_start,width)
	return total_area

def draw_trapezium(ax,x0,x1,y0,y1,y2,mycolor='blue'):
	verts = [
		(x0, y0), # left, bottom
		(x0, y1), # left, top
		(x1, y2), # right, top
		(x1, y0), # right, bottom
		(0., 0.), # ignored
	]

	codes = [Path.MOVETO,
		 Path.LINETO,
		 Path.LINETO,
		 Path.LINETO,
		 Path.CLOSEPOLY,
		 ]

	path = Path(verts, codes)

	patch = patches.PathPatch(path, facecolor='none',alpha=0.6,edgecolor=mycolor,ls='solid',lw=2)
	ax.add_patch(patch)

def draw_rectangle(ax,x0,x1,y0,y1,mycolor='green'):
	draw_trapezium(ax,x0,x1,y0,y1,y1,mycolor)

def draw_trapezium_strips(t_start, t_end, n_strips,ax):
	interval=t_end-t_start
	width=interval/n_strips
	total_area=0
	for i in range(n_strips):
		strip_start=t_start+ i *width
		strip_end=strip_start+width
		draw_trapezium(ax,strip_start,strip_end,0,dydt(strip_start),dydt(strip_end))
	return True

def draw_rectangular_strips(t_start, t_end, n_strips,ax):
	interval=t_end-t_start
	width=interval/n_strips
	total_area=0
	for i in range(n_strips):
		strip_start=t_start+ i *width
		strip_end=strip_start+width
		draw_rectangle(ax,strip_start,strip_end,0,dydt(strip_start))
	return True


#########################################
#
# You can change the parameters used
# for integration below.
#

n_strips = 5
t_start = 0.0
t_end = 3.0

# If you need to change the y-axis limits change the lines below
#
y_axis_min = 0
y_axis_max = 20


##########################################
# DO NOT CHANGE THE CODE BELOW HERE....
##########################################

print ("The exact value of the integral is:")
print (y(t_end) - y(t_start))

print ("Using {} rectangles to approximate the area we calculate the integral as:".format(n_strips))
print (sum_strips(t_start,t_end,n_strips))

#print "Using {} trapezoids to approximate the area we calculate the integral as:".format(n_strips)
#print sum_trapezium_strips(t_start,t_end,n_strips)

# COD TO DRAW FIGURE SHOWING STRIPS METHOD
plt.close("all")
fig=plt.figure(figsize=[10,6])

t_full=np.linspace(t_start-1.,t_end+1.,1000)
y_full=[]

dydt_full=[]
t_shaded=[]
dydt_shaded=[]

for ti in t_full:
	dydt_i=dydt(ti)
	dydt_full.append(dydt_i)
	y_full.append(y(ti))
	if t_start <= ti <= t_end:
		t_shaded.append(ti)
		dydt_shaded.append(dydt_i)

ax=fig.add_subplot(1,1,1)
ax.plot(t_full,dydt_full,ls='dotted')
ax.fill_between(t_shaded,dydt_shaded,alpha=0.2)
ax.text(-0.5,10,"Numerical intergration to \nestimate area under curve of \ny'(t) = t^3 - 3*t^2/2 + 4")

draw_rectangular_strips(t_start,t_end,n_strips,ax)
#draw_trapezium_strips(t_start,t_end,n_strips,ax)

ax.set_ylim(0,y_axis_max)
ax.set_ylabel("dy/dt")

fig.show()
