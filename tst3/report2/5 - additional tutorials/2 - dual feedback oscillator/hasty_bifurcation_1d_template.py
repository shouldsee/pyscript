import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

def sdot(s,t,p):
    (v,K,T,kd_A,kd_R)=p
    A,R=s
    dA =
    dR =
    return (dA,dR)

# set parameters:
k_dA=
k_dR=

v=
K=
Z=

# set initial conditions
A0=
R0=

s0=(A0,R0)
params=(v,K,Z,k_dA,k_dR)

t_max=2000
n_obs=t_max*10
t_obs=np.linspace(0,t_max,n_obs)

plt.close("all")
fig1 = plt.figure()
ax=plt.subplot(1,1,1)

k_dR_list=np.linspace(0.01,0.08,8)
for k_dR in k_dR_list:
    params=(v,K,Z,k_dA,k_dR)
    s_obs=odeint(sdot,s0,t_obs,args=(params,))

    # we will slice output to only include observations
    # from second half of the simulation
    # so transient behaviour is not included
    midpoint=int(n_obs/2)

    # extract observations of A and R from midpoint onwards
    A_obs=s_obs[midpoint:,0]
    R_obs=s_obs[midpoint:,1]

    max_A=max(A_obs)
    min_A=
    range_A=

    mystr="Input: k_dA={}, k_dR={:.3f}.  Output [A]: max={:5.2f}, min={:5.2f}, range={:.5f}"
    print mystr.format(k_dA,k_dR, max_A,min_A,range_A)

fig1.show()
