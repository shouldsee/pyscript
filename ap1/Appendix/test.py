import numpy as np
import math as m
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from functools import reduce
import csv
#%matplotlib inline
#from utils import *
import copy

def fNdot(N,t,params):
    (s,i,r)=N;
    T,r0=params;
    ds=0.      -r0*s*i/T;
    di=r0*s*i/T-i/T;
    dr=i/T;
    return (ds,di,dr)

def vfield(axis,Ndot,num=20,scale=1/0.0015,norm_method='none',col=None):
    xlim=list(axis.get_xlim());
    ylim=list(axis.get_ylim());
    x=np.linspace(xlim[0],xlim[1],num);
    y=np.linspace(ylim[0],ylim[1],num);
    xs,ys=np.meshgrid(x,y);
    us,vs=Ndot(xs,ys);
    
    ls=(us**2+vs**2)**0.5;        
    if norm_method=='log':
        lls=np.min(np.log(ls),0);
        us=us/ls*lls;
        vs=vs/ls*lls;
    elif norm_method=='len':
        ls=(us**2+vs**2)**0.5;
        us=us/ls;
        vs=vs/ls;
    axis.quiver(xs,ys,us,vs,color=col,scale_units='inches',angles='xy',pivot='mid',scale=scale)



class intobj():
    def __init__(self,Ndot,N0,params):
        self.fcn=Ndot;
        self.ss,self.labels=N0;
        self.params=params;
    def evo(self,ts):
        self.ts=ts;
        Ns=odeint(self.fcn,self.ss,ts,args=(self.params,));
        self.Ns=Ns;
        return(Ns)
    def line(self,ax,idx):
        ax.plot(self.ts,self.Ns[:,idx],label=self.labels[idx]);
    def phase(self,ax,idx):
        xi,yi=idx;
        ax.plot(self.Ns[:,xi],self.Ns[:,yi]);
        ax.set_xlabel(self.labels[xi]);
        ax.set_ylabel(self.labels[yi]);
        


from scipy.interpolate import *

def readcsv(fname):
#     fname='measles_data.csv'
    csvdata={};

    with open(fname,'r') as  f:
        header=f.readline().rstrip('\n').split(',');
        lsts=[[]]*(len(header));
        for line in f.readlines():
            data=line.rstrip('\n').split(',');
            for i,k in enumerate(data):
    #             len(lsts)
                lsts[i]=lsts[i]+[float(k)];
    #     print(lsts)
        for i,lst in enumerate(lsts):
            csvdata[header[i]]=lst;
    #     print(len(lsts))
    return csvdata;
csvdata=readcsv('measles_data.csv')
time_years=csvdata['time'];
time_days=[round(x*365) for x in csvdata['time']];
Is=csvdata['infection'];


fname='EWpop_short.csv'
popcsv=readcsv(fname);
pop_pre=popcsv['pop'];
f=interp1d(list(x + 0.5 for x in popcsv['year']),pop_pre);
pop=f(time_years);


### figure 1
### visualise historical data
time_years=csvdata['time'];
time_days=[round(x*365) for x in csvdata['time']];
infection=csvdata['infection'];
Is=[x/y for x,y in zip(infection,pop)];
fig=plt.figure(figsize=[10,8]);
ax=plt.subplot(2,1,1);
ax.plot(time_years,Is);
ax.set_title ('Fraction of England and Wales population infected with Measles over 10 years')
ax.set_xlabel('time/year')
ax.set_ylabel('frac_I')

# ax2.plot(time_years[0:1500],Is[0:1500]);
# ax2.set_title ('Fraction of England and Wales population infected with Measles over 1500 days')
# ax2.set_xlabel('time/year')
# ax2.set_ylabel('frac_I')

fig.show()     