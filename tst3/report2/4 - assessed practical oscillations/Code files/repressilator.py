import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
import math

def sdot(s,t,param):

    (p_LacI, p_TetR, p_CI, m_LacI, m_TetR, m_CI, m_GFP, p_GFP)=s
    (k_m, k_m0, k_p, k_dm, k_dp, K, n, k_dGFP)=param

    rate_m_LacI_prod = k_m*K**n / (K**n + p_CI**n)   + k_m0
    rate_m_TetR_prod = k_m*K**n / (K**n + p_LacI**n) + k_m0
    rate_m_CI_prod   = k_m*K**n / (K**n + p_TetR**n) + k_m0
    rate_m_GFP_prod  = 0

    rate_p_LacI_prod = k_p*m_LacI
    rate_p_TetR_prod = k_p*m_TetR
    rate_p_CI_prod   = k_p*m_CI
    rate_p_GFP_prod  = 0

    rate_m_LacI_loss = k_dm*m_LacI
    rate_m_TetR_loss = k_dm*m_TetR
    rate_m_CI_loss   = k_dm*m_CI
    rate_m_GFP_loss  = 0

    rate_p_LacI_loss = k_dp*p_LacI
    rate_p_TetR_loss = k_dp*p_TetR
    rate_p_CI_loss   = k_dp*p_CI
    rate_p_GFP_loss  = 0

    dp_LacI = rate_p_LacI_prod - rate_p_LacI_loss
    dp_TetR = rate_p_TetR_prod - rate_p_TetR_loss
    dp_CI   = rate_p_CI_prod   - rate_p_CI_loss
    dp_GFP  = rate_p_GFP_prod  - rate_p_GFP_loss


    dm_LacI = rate_m_LacI_prod - rate_m_LacI_loss
    dm_TetR = rate_m_TetR_prod - rate_m_TetR_loss
    dm_CI   = rate_m_CI_prod   - rate_m_CI_loss
    dm_GFP  = rate_m_GFP_prod   - rate_m_GFP_loss

    sdot = (dp_LacI, dp_TetR, dp_CI, dm_LacI, dm_TetR, dm_CI, dm_GFP, dp_GFP)
    return sdot

# DEFINE INITIAL CONDITIONS AND PARAMETERS

#intitial condtions
p_LacI=0
p_TetR=0
p_CI=0
p_GFP=0

m_LacI=5
m_TetR=0
m_CI=0
m_GFP=0

s = (p_LacI, p_TetR, p_CI, m_LacI, m_TetR, m_CI)

k_m=0.5
k_m0=5e-4

K=40.0
n=2.1

k_dm=0.00577622650467
k_dp=0.00144405662617
k_dGFP=0.000128360588993
k_p=0.115524530093

s0 = (p_LacI, p_TetR, p_CI, m_LacI, m_TetR, m_CI, m_GFP, p_GFP)
param = (k_m, k_m0, k_p, k_dm, k_dp, K, n, k_dGFP)

t_max=1000.*60.
t_obs=np.linspace(0,t_max,t_max+1)
s_obs=odeint(sdot,s0,t_obs,args=(param,))

p_LacI_obs = s_obs[:,0]
p_TetR_obs = s_obs[:,1]
p_CI_obs   = s_obs[:,2]
m_LacI_obs = s_obs[:,3]
m_TetR_obs = s_obs[:,4]
m_CI_obs   = s_obs[:,5]
m_GFP_obs   = s_obs[:,6]
p_GFP_obs   = s_obs[:,7]
