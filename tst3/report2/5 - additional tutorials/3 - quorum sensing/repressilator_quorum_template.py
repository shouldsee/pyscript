import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
import math

def sdot(s,t,param):

    (k_m, k_m0, k_p, k_dm, k_dp, K, n, K, K2, k_m2, k_diff, V, k_dAI_ext )=param

    (p_LacI_1, p_TetR_1, p_CI_1, p_AI_1, m_LacI_1, m_TetR_1, m_CI_1, m_AI_1,
     p_LacI_2, p_TetR_2, p_CI_2, p_AI_2, m_LacI_2, m_TetR_2, m_CI_2, m_AI_2,
     p_AI_ext) = s

    # CELL 1
    rate_m_LacI_prod_1 = k_m*K**n / (K**n + p_CI_1**n)   + k_m2 * p_AI_1 / (K2 + p_AI_1) + k_m0
    rate_m_TetR_prod_1 = k_m*K**n / (K**n + p_LacI_1**n) + k_m0
    rate_m_CI_prod_1   = k_m*K**n / (K**n + p_TetR_1**n) + k_m0
    rate_m_AI_prod_1   = k_m*K**n / (K**n + p_LacI_1**n) + k_m0

    rate_p_LacI_prod_1 = k_p*m_LacI_1
    rate_p_TetR_prod_1 = k_p*m_TetR_1
    rate_p_CI_prod_1   = k_p*m_CI_1
    rate_p_AI_prod_1   = k_p*m_AI_1

    rate_m_LacI_loss_1 = k_dm*m_LacI_1
    rate_m_TetR_loss_1 = k_dm*m_TetR_1
    rate_m_CI_loss_1   = k_dm*m_CI_1
    rate_m_AI_loss_1   = k_dm*m_AI_1

    rate_p_LacI_loss_1 = k_dp*p_LacI_1
    rate_p_TetR_loss_1 = k_dp*p_TetR_1
    rate_p_CI_loss_1   = k_dp*p_CI_1
    rate_p_AI_loss_1   = k_dp*p_AI_1

    rate_p_AI_diff_net_1  = 0

    dp_LacI_1 = rate_p_LacI_prod_1 - rate_p_LacI_loss_1
    dp_TetR_1 = rate_p_TetR_prod_1 - rate_p_TetR_loss_1
    dp_CI_1   = rate_p_CI_prod_1   - rate_p_CI_loss_1
    dp_AI_1   = rate_p_AI_prod_1   - rate_p_AI_loss_1

    dm_LacI_1 = rate_m_LacI_prod_1 - rate_m_LacI_loss_1
    dm_TetR_1 = rate_m_TetR_prod_1 - rate_m_TetR_loss_1
    dm_CI_1   = rate_m_CI_prod_1   - rate_m_CI_loss_1
    dm_AI_1   = rate_m_AI_prod_1   - rate_m_AI_loss_1


    # CELL 2
    rate_m_LacI_prod_2 = k_m*K**n / (K**n + p_CI_2**n)   + k_m2 * p_AI_2 / (K2 + p_AI_2) + k_m0
    rate_m_TetR_prod_2 = k_m*K**n / (K**n + p_LacI_2**n) + k_m0
    rate_m_CI_prod_2   = k_m*K**n / (K**n + p_TetR_2**n) + k_m0
    rate_m_AI_prod_2   = k_m*K**n / (K**n + p_LacI_2**n) + k_m0

    rate_p_LacI_prod_2 = k_p*m_LacI_2
    rate_p_TetR_prod_2 = k_p*m_TetR_2
    rate_p_CI_prod_2   = k_p*m_CI_2
    rate_p_AI_prod_2   = k_p*m_AI_2

    rate_m_LacI_loss_2 = k_dm*m_LacI_2
    rate_m_TetR_loss_2 = k_dm*m_TetR_2
    rate_m_CI_loss_2   = k_dm*m_CI_2
    rate_m_AI_loss_2   = k_dm*m_AI_2

    rate_p_LacI_loss_2 = k_dp*p_LacI_2
    rate_p_TetR_loss_2 = k_dp*p_TetR_2
    rate_p_CI_loss_2   = k_dp*p_CI_2
    rate_p_AI_loss_2   = k_dp*p_AI_2

    rate_p_AI_diff_net_2  = 0

    dp_LacI_2 = rate_p_LacI_prod_2 - rate_p_LacI_loss_2
    dp_TetR_2 = rate_p_TetR_prod_2 - rate_p_TetR_loss_2
    dp_CI_2   = rate_p_CI_prod_2   - rate_p_CI_loss_2
    dp_AI_2   = rate_p_AI_prod_2   - rate_p_AI_loss_2

    dm_LacI_2 = rate_m_LacI_prod_2 - rate_m_LacI_loss_2
    dm_TetR_2 = rate_m_TetR_prod_2 - rate_m_TetR_loss_2
    dm_CI_2   = rate_m_CI_prod_2   - rate_m_CI_loss_2
    dm_AI_2   = rate_m_AI_prod_2   - rate_m_AI_loss_2

    rate_p_AI_ext_loss = k_dAI_ext * p_AI_ext

    dp_AI_ext =  - rate_p_AI_ext_loss

    ds = (dp_LacI_1, dp_TetR_1, dp_CI_1, dp_AI_1, dm_LacI_1, dm_TetR_1, dm_CI_1, dm_AI_1,
        dp_LacI_2, dp_TetR_2, dp_CI_2, dp_AI_2, dm_LacI_2, dm_TetR_2, dm_CI_2, dm_AI_2,
    dp_AI_ext)
    return ds

# DEFINE INITIAL CONDITIONS AND PARAMETERS

#intitial condtions
p_LacI_1=0
p_TetR_1=0
p_CI_1=0
p_AI_1=0

m_LacI_1=4
m_TetR_1=0
m_CI_1=0
m_AI_1=0

p_LacI_2=0
p_TetR_2=0
p_CI_2=0
p_AI_2=0

m_LacI_2=0
m_TetR_2=4
m_CI_2=0
m_AI_2=0
V=5
p_AI_ext=1.0

# paramters

k_m=0.5
k_m0=5e-4
t_half_p=8.*60
t_half_m=2.*60
av_p_per_mRNA=20.0
K=40.0
n=2.1

k_dm=math.log(2.0)/t_half_m
k_dp=math.log(2.0)/t_half_p

t_av_p=1./k_dp
t_av_m=1./k_dm
k_p=av_p_per_mRNA/t_av_m

k_diff=math.log(2.0)/(2*60)
k_dAI_ext=math.log(2.0)/(30*60)
k_m2=0.01
K2=100.0

s0=(p_LacI_1, p_TetR_1, p_CI_1, p_AI_1, m_LacI_1, m_TetR_1, m_CI_1, m_AI_1,
    p_LacI_2, p_TetR_2, p_CI_2, p_AI_2, m_LacI_2, m_TetR_2, m_CI_2, m_AI_2,
    p_AI_ext)
param=    (k_m, k_m0, k_p, k_dm, k_dp, K, n, K, K2, k_m2, k_diff, V, k_dAI_ext)

t_max=1500.*60.
t_obs=np.linspace(0,t_max,t_max+1)
s_obs=odeint(sdot,s0,t_obs,args=(param,))
plt.close('all')

p_LacI_1_obs = s_obs[:,0]
p_TetR_1_obs = s_obs[:,1]
p_CI_1_obs   = s_obs[:,2]
p_AI_1_obs   = s_obs[:,3]
m_LacI_1_obs = s_obs[:,4]
m_TetR_1_obs = s_obs[:,5]
m_CI_1_obs   = s_obs[:,6]
m_AI_1_obs   = s_obs[:,7]

p_LacI_2_obs = s_obs[:,8]
p_TetR_2_obs = s_obs[:,9]
p_CI_2_obs   = s_obs[:,10]
p_AI_2_obs   = s_obs[:,11]
m_LacI_2_obs = s_obs[:,12]
m_TetR_2_obs = s_obs[:,13]
m_CI_2_obs   = s_obs[:,14]
m_AI_2_obs   = s_obs[:,15]

p_AI_ext = s_obs[:,16]

fig=plt.figure()
# add "super title" over plot set
fig.suptitle('Simple Gene Expression Model')
ax_P=fig.add_subplot(3,1,1)
ax_mRNA=fig.add_subplot(3,1,2)
ax_AI=fig.add_subplot(3,1,3)

ax_P.plot(t_obs/60, p_LacI_1_obs, 'b-')
ax_P.plot(t_obs/60, p_LacI_2_obs, 'r-')

ax_mRNA.plot(t_obs/60, m_LacI_1_obs, 'b-')
ax_mRNA.plot(t_obs/60, m_LacI_2_obs, 'r-')

ax_AI.plot(t_obs/60, p_AI_1_obs, 'b-')
ax_AI.plot(t_obs/60, p_AI_2_obs, 'r-')
ax_AI.plot(t_obs/60, p_AI_ext, 'k-')

ax_P.set_ylabel('protein LacI')
ax_mRNA.set_ylabel('mRNA LacI')
ax_AI.set_ylabel('AI level')

ax_AI.set_xlabel('Time / min')
fig.show()
