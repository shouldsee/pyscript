import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

# MODEL FOR DETERMINISTIC GENE EXPRESSION

# SPECIES IN STATE VECTOR

# s=( g_on,  # conc active gene
#     g_off, # conc inactive gene
#     mRNA,  # conc RNA
#     p      # conc protien
#    )


# PARAMETERS USED

# param=( k_on    #rate_constant_on,
#         k_off,  #rate_constant_off,
#         k_t,    #rate_constant_transcription,
#         k_p,    # rate_constant_translation
#         kd_mRNA # degradation rate constant mRNA
#         kd_p    # degradation/dilution rate constant protein
#         )


def sdot_gene_expression(species,t,param):

    g_on,g_off,mRNA,p = species
    k_on, k_off, k_t, k_p, kd_mRNA, kd_p = param

    rate_on = k_on*g_off
    rate_off =
    rate_t =
    rate_p =

    rate_dmRNA =
    rate_dp =

    dg_on = rate_on - rate_off
    dg_off =
    dmRNA =
    dp =

    sdot = (dg_on,dg_off,dmRNA,dp)
    return sdot


# DEFINE INITIAL CONDITIONS AND PARAMETERS

g_off=
g_on=
mRNA=
p=

k_on=
k_off=
k_t=
k_p=
kd_mRNA=
kd_p=

s0=(g_on,g_off,mRNA,p)
param=(k_on, k_off, k_t, k_p, kd_mRNA, kd_p)

t_max=
t_obs=np.linspace(0,t_max,10001)
s_obs=odeint(sdot_gene_expression,s0,t_obs,args=(param,))

fig=plt.figure()
# add "super title" over plot set
fig.suptitle('Simple Gene Expression Model')
ax_G=fig.add_subplot(3,1,1)
ax_mRNA=fig.add_subplot(3,1,2)
ax_P=fig.add_subplot(3,1,3)

ax_G.set_ylabel('fraction $g_{on}$')
ax_mRNA.set_ylabel('#mRNA copies')
ax_P.set_ylabel('#Protein copies')
ax_P.set_xlabel('time (s)')

ax_G.plot(t_obs, s_obs[:,0], 'b-')
ax_mRNA.plot(t_obs, s_obs[:,2], 'r-')
ax_P.plot(t_obs, s_obs[:,3], 'g-')

# adjusts matplotloib to autoscale margins with
# 0.1 i.e 10% head room at the top of the y-axis
ax_P.margins(y=0.1)
ax_P.set_ylim(bottom=0)
ax_mRNA.margins(y=0.1)
ax_mRNA.set_ylim(bottom=0)

fig.show()
