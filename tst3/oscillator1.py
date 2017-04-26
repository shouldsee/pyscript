import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
import math

def sdot_oscillator(s,t,p):
    (X,V) = s
    (k1,k2,k3) = p
    dX = V
    dV = -k1 * (X - 8.)
    return (dX,dV)

# set parameter values
# DO NOT EDIT!
T=5.0
k1=(2.0*math.pi/T)**2
k2=0.2
k3=0.1
p=(k1,k2,k3)

# set up figure
plt.close("all")
fig1 = plt.figure(figsize=(14,6))
axL=fig1.add_subplot(1,2,1)
axL.set_xlabel('t')
axL.set_ylabel('X')
axL.set_xlim(0,35)
axL.set_ylim(0,20)
axR=fig1.add_subplot(1,2,2)


# simulate oscillator
t_obs=np.linspace(0,35,10001)
s0=[4.0, 0.0]
s_obs=odeint(sdot_oscillator,s0,t_obs,args=(p,))


# show figure


lbl='traj1';
axL.plot(t_obs,s_obs[:,0],label=lbl)
axL.set_xlabel('t')
axL.set_ylabel('X')

axR.plot(s_obs[:,0],s_obs[:,1],label=lbl);
axR.set_xlabel('X')
axR.set_ylabel('Y')

s0=[8.3, 0.0]
s_obs=odeint(sdot_oscillator,s0,t_obs,args=(p,))
lbl='traj2';
axL.plot(t_obs,s_obs[:,0],label=lbl)
axL.set_xlabel('t')
axL.set_ylabel('X')

axR.plot(s_obs[:,0],s_obs[:,1],label=lbl);
axR.set_xlabel('X')
axR.set_ylabel('Y')

axL.legend();
axR.legend();
fig1.show()

fig1.savefig('oscillator1.png')
