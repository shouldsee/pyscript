import numpy as np
from matplotlib import pyplot as plt
from scipy.integrate import odeint
from matplotlib.widgets import Slider, Button, RadioButtons

def sdot(s,t,params):
    k1,KI,n,k2,k3,k4,k5,k6,k7,k8=params
    X,Y1,Y2,Z=s
    dX = k1/(KI**n + Z**n) - k2*X
    dY1 = k3*X - k4*Y1
    dY2 = k7*Y1 - k8*Y2
    dZ = k5*Y2 - k6*Z
    ds=(dX,dY1,dY2,dZ)
    return ds

k1=1.0
KI=1.0
n=5.0
k2=1.0
k3=10.0
k4=1.0
k5=10.0
k6=1.0
k7=10.0
k8=1.0
params=k1,KI,n,k2,k3,k4,k5,k6,k7,k8

X=0.0
Y1=0.0
Y2=0.0
Z=0.00
s0=X,Y1,Y2,Z

## RUN MODEL WITH DEFAULT PARAMETER VALUES
t_max=1000
t_obs=np.linspace(0,t_max,t_max*10)
s_obs=odeint(sdot,s0,t_obs,args=(params,))

# trim transient behaviour
t_trim=800
obs_trim=int(t_trim*10) # because we used linspace with t_max*10 intervals
t_obs=t_obs[obs_trim:]
s_obs=s_obs[obs_trim:]

## CREATE TIMESERIES PLOT
plt.close("all")
fig_timeseries = plt.figure(figsize=(12,8))
ax1=fig_timeseries.add_subplot(3,1,1)
ax2=fig_timeseries.add_subplot(3,1,2)
ax3=fig_timeseries.add_subplot(3,1,3)

fig_timeseries.suptitle('Concentration vs Time')
ax1.plot(t_obs, s_obs[:,0], '-',label='[X]')
ax2.plot(t_obs, s_obs[:,1], '-',label='[Y]')
ax3.plot(t_obs, s_obs[:,3], '-',label='[Z]')
ax1.legend()
ax2.legend()
ax3.legend()
fig_timeseries.show()

## CREATE PHASE PLOT Z vs X
fig_phase = plt.figure(figsize=(12,8))
ax=fig_phase.add_subplot(1,1,1)
ax.plot(s_obs[:,0],s_obs[:,3])
ax.set_title('Phase space Z vs X')
ax.set_xlabel('X')
ax.set_ylabel('Z')

fig_phase.show()
