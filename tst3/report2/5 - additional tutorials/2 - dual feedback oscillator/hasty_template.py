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

t_max=800
n_obs=t_max*10
t_obs=np.linspace(0,t_max,n_obs)
s_obs=odeint(sdot,s0,t_obs,args=(params,))
A_obs=s_obs[:,0]
R_obs=s_obs[:,1]
