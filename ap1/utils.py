import numpy as np
import math as m
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from functools import reduce
import csv
def vfield(axis,Ndot,num=20,scale=1/0.0015,norm_method='none'):
    xlim=list(axis.get_xlim());
    ylim=list(axis.get_ylim());
    x=np.linspace(xlim[0],xlim[1],num);
    y=np.linspace(ylim[0],ylim[1],num);
    xs,ys=np.meshgrid(x,y);
    us,vs=Ndot((xs,ys),0,params);
    
    ls=(us**2+vs**2)**0.5;        
    if norm_method=='log':
        lls=np.min(np.log(ls),0);
        us=us/ls*lls;
        vs=vs/ls*lls;
    elif norm_method=='len':
        ls=(us**2+vs**2)**0.5;
        us=us/ls;
        vs=vs/ls;
    axis.quiver(xs,ys,us,vs,color='b',scale_units='inches',angles='xy',pivot='mid',scale=scale)

class intobj():
    def __init__(self,Ndot,N0,params):
        self.fcn=Ndot;
        self.ss,self.labels=N0;
        self.params=params;
    def evo(self,ts):
        self.ts=ts;
        Ns=odeint(self.fcn,self.ss,ts,args=(self.params,));
        intobj.Ns=Ns;
        return(Ns)
    def line(self,ax,idx,label=[]):
	if label==[]:
		label=self.labels[idx]
	else:
		pass
        ax.plot(self.ts,self.Ns[:,idx],label=label);
    def phase(self,ax,idx):
        xi,yi=idx;
        xlabel=self.labels[xi];
        ax.plot(self.Ns[:,xi],self.Ns[:,yi]);
        # ax.set_xlabel(self.labels[xi]);
        # ax.set_ylabel(self.labels[yi]);
        # ax.title('Phase plot for ',self.labels[xi],' and ',self.labels[yi])


def fNdot(N,t,params):
#     tr=T, r0=r0;
    (s,i,r)=N;
    T,r0=params;
    ds=0.      -r0*s*i/T;
    di=r0*s*i/T-i/T;
    dr=i/T;
    return (ds,di,dr)
