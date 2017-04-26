import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
import math

def sdot_oscillator(s,t,p):
    (X,V) = s
    (k1,k2,k3) = p
    dX =
    dV =
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
fig1.show()
