import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
import math

# function that models a transcription rate
# dependant on the night/day light level
# with a 24 hour cycle
def L(t):
    return ( math.sin(math.pi*t/24) )**2.0

def sdot(s,t,params):

    M,PC,PN=s
    vs,vm,Km,ks,vd,k1,k2,KI,Kd,n,xmax=params

    dM =
    dPC =
    dPN =

    return dM,dPC,dPN

# initialise parameter values

vs=
vm=
Km=
ks=
vd=
k1=
k2=
KI=
Kd=
n=

# set initial conditions

M=
P0=
PN=

xmax=0.0

s0=M,P0,PN
params=vs,vm,Km,ks,vd,k1,k2,KI,Kd,n,xmax

# simulate system for 12 days
t_max=24*12
t_obs=np.linspace(0,t_max,t_max*10)
s_obs=odeint(sdot,s0,t_obs,args=(params,))

# create some output plots
plt.close("all")

# show plot of mRNA and protein levels
fig1 = plt.figure()
ax=plt.subplot(1,1,1)
ax.set_title('Simple model of a circadian clock gene')
ax.plot(t_obs/24,s_obs[:,0],label='mRNA') # we divide t by 24 to convert hours to days
ax.plot(t_obs/24,s_obs[:,1],label='protein (cellular)')
ax.plot(t_obs/24,s_obs[:,2],label='protein (nuclear)')
ax.set_ylabel('conc. (arb)')
ax.set_xlabel('time / days')
fig1.show()


# generate a plot of the 24 hour light cycle
L_obs=[]
for t in t_obs:
    L_obs.append(L(t))

# compare mRNA level to 24 hour light cycle
fig2 = plt.figure()
fig2.suptitle('Simple model of a circadian clock gene (coupling parameter xmax={})'.format(xmax))
# top plot shows mRNA level
axH=plt.subplot(2,1,1)
axH.set_ylabel('mRNA conc. (arb)')
axH.plot(t_obs/24,s_obs[:,0],'r') # divide t by 24 to convert hours to days
# lower plot shows light cycle
axL=plt.subplot(2,1,2)
axL.plot(t_obs/24,L_obs,'b') # divide t by 24 to convert hours to days
axL.set_ylabel('Light level (arb)')
axL.set_xlabel('time / days')
fig2.show()
