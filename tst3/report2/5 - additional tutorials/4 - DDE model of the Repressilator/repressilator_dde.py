import numpy as np
import scipy.integrate
from scipy.integrate import odeint
import scipy.interpolate
from matplotlib import pyplot as plt
import math

class ddeVar:
    """ special function-like variables for the integration of DDEs """


    def __init__(self,g,tc=0):
        """ g(t) = expression of Y(t) for t<tc """

        self.g = g
        self.tc= tc
        # We must fill the interpolator with 2 points minimum
        self.itpr = scipy.interpolate.interp1d(
            np.array([tc-1,tc]), # X
            np.array([self.g(tc),self.g(tc)]).T, # Y
            kind='linear', bounds_error=False,
            fill_value = self.g(tc))


    def update(self,t,Y):
        """ Add one new (ti,yi) to the interpolator """

        self.itpr.x = np.hstack([self.itpr.x, [t]])
        Y2 = Y if (Y.size==1) else np.array([Y]).T
        self.itpr.y = np.hstack([self.itpr.y, Y2])
        self.itpr._y = self.itpr._reshape_yi(self.itpr.y)
        self.itpr.fill_value = Y


    def __call__(self,t=0):
        """ Y(t) will return the instance's value at time t """
        #print "t",t
        #print self.itpr(t)
        return (self.g(t) if (t<=self.tc) else self.itpr(t))



class dde(scipy.integrate.ode):
    """ Overwrites a few functions of scipy.integrate.ode"""


    def __init__(self,f,jac=None):

        def f2(t,y,args):
            return f(self.Y,t,*args)
        scipy.integrate.ode.__init__(self,f2,jac)
        self.set_f_params(None)


    def integrate(self, t, step=0, relax=0):

        scipy.integrate.ode.integrate(self,t,step,relax)
        self.Y.update(self.t,self.y)
        return self.y


    def set_initial_value(self,Y):

        self.Y = Y #!!! Y will be modified during integration
        scipy.integrate.ode.set_initial_value(self, Y(Y.tc), Y.tc)



def ddeint(func,g0,tt,args=None):
    """ similar to scipy.integrate.odeint. Solves the DDE system
        defined by func at the times tt with 'history function' g
        and potential additional arguments for the model, fargs
    """
    dde_ = dde(func)
    g=lambda t:  np.array(g0)
    dde_.set_initial_value(ddeVar(g,tt[0]))
    dde_.set_f_params(args if args else [])
    return np.array([g(tt[0])]+[dde_.integrate(dde_.t + dt)
                                 for dt in np.diff(tt)])

# rate function for ODE modelling of Repressilator
def sdot_ODE(s,t,param):

    # load levels of proteins/mRNA
    (p_LacI, p_TetR, p_CI, m_LacI, m_TetR, m_CI) = s

    # load parameter values
    (k_m, k_m0, k_p, k_dm, k_dp, K, n )=param

    rate_m_LacI_prod = k_m*K**n / (K**n + p_CI**n)   + k_m0
    rate_m_TetR_prod = k_m*K**n / (K**n + p_LacI**n) + k_m0
    rate_m_CI_prod   = k_m*K**n / (K**n + p_TetR**n) + k_m0

    rate_p_LacI_prod = k_p*m_LacI
    rate_p_TetR_prod = k_p*m_TetR
    rate_p_CI_prod   = k_p*m_CI

    rate_m_LacI_loss = k_dm*m_LacI
    rate_m_TetR_loss = k_dm*m_TetR
    rate_m_CI_loss   = k_dm*m_CI

    rate_p_LacI_loss = k_dp*p_LacI
    rate_p_TetR_loss = k_dp*p_TetR
    rate_p_CI_loss   = k_dp*p_CI

    dp_LacI = rate_p_LacI_prod - rate_p_LacI_loss
    dp_TetR = rate_p_TetR_prod - rate_p_TetR_loss
    dp_CI   = rate_p_CI_prod   - rate_p_CI_loss

    dm_LacI = rate_m_LacI_prod - rate_m_LacI_loss
    dm_TetR = rate_m_TetR_prod - rate_m_TetR_loss
    dm_CI   = rate_m_CI_prod   - rate_m_CI_loss

    ds = [dp_LacI, dp_TetR, dp_CI, dm_LacI, dm_TetR, dm_CI]
    return ds

# rate function for DDE modelling of Repressilator
def sdot_DDE(s,t,param):

    # load current levels of proteins/mRNA
    (p_LacI, p_TetR, p_CI, m_LacI, m_TetR, m_CI) = s(t)

    # transcription / translation delay both set as 60 seconds
    delay = 60

    # load past levels of proteins/mRNA present in system 60 seconds previously
    (pd_LacI, pd_TetR, pd_CI, md_LacI, md_TetR, md_CI) = s(t-delay)

    # load parameter values
    (k_m, k_m0, k_p, k_dm, k_dp, K, n )=param

    # lines 136 to 138 and 140 to 142 need editing so that
    # the past values of protein and mRNA levels are used
    # where appropriate (see tutorial for guidance)

    rate_m_LacI_prod = k_m*K**n / (K**n + p_CI**n)   + k_m0
    rate_m_TetR_prod = k_m*K**n / (K**n + p_LacI**n) + k_m0
    rate_m_CI_prod   = k_m*K**n / (K**n + p_TetR**n) + k_m0

    rate_p_LacI_prod = k_p*m_LacI
    rate_p_TetR_prod = k_p*m_TetR
    rate_p_CI_prod   = k_p*m_CI

    rate_m_LacI_loss = k_dm*m_LacI
    rate_m_TetR_loss = k_dm*m_TetR
    rate_m_CI_loss   = k_dm*m_CI

    rate_p_LacI_loss = k_dp*p_LacI
    rate_p_TetR_loss = k_dp*p_TetR
    rate_p_CI_loss   = k_dp*p_CI

    dp_LacI = rate_p_LacI_prod - rate_p_LacI_loss
    dp_TetR = rate_p_TetR_prod - rate_p_TetR_loss
    dp_CI   = rate_p_CI_prod   - rate_p_CI_loss

    dm_LacI = rate_m_LacI_prod - rate_m_LacI_loss
    dm_TetR = rate_m_TetR_prod - rate_m_TetR_loss
    dm_CI   = rate_m_CI_prod   - rate_m_CI_loss

    sdot = [dp_LacI, dp_TetR, dp_CI, dm_LacI, dm_TetR, dm_CI]
    return sdot

# DEFINE INITIAL CONDITIONS AND PARAMETERS

#intitial condtions
p_LacI=0
p_TetR=0
p_CI=0

m_LacI=5
m_TetR=0
m_CI=0

# parameters
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


s0=(p_LacI, p_TetR, p_CI, m_LacI, m_TetR, m_CI)
param=    (k_m, k_m0, k_p, k_dm, k_dp, K, n)

t_max=600.*60.
t_obs=np.linspace(0,t_max,t_max+1)

print "Running ODE"
s_obs1=odeint(sdot_ODE,s0,t_obs,args=(param,))
print "Finished ODE"

print "Running DDE"
s_obs2=ddeint(sdot_DDE,s0,t_obs,args=(param,))
print "Finished DDE"

plt.close('all')

p_LacI_obs1 = s_obs1[:,0]
p_TetR_obs1 = s_obs1[:,1]
p_CI_obs1   = s_obs1[:,2]

p_LacI_obs2 = s_obs2[:,0]
p_TetR_obs2 = s_obs2[:,1]
p_CI_obs2   = s_obs2[:,2]

fig=plt.figure()
# add "super title" over plot set
fig.suptitle('Repressilator Model (ODE top vs DDE bottom)')
ax_ODE=fig.add_subplot(2,1,1)
ax_DDE=fig.add_subplot(2,1,2)

ax_ODE.plot(t_obs/60, p_LacI_obs1, 'b-')
ax_ODE.plot(t_obs/60, p_TetR_obs1, 'r-')
ax_ODE.plot(t_obs/60, p_CI_obs1, 'g-')

ax_DDE.plot(t_obs/60, p_LacI_obs2, 'b-')
ax_DDE.plot(t_obs/60, p_TetR_obs2, 'r-')
ax_DDE.plot(t_obs/60, p_CI_obs2, 'g-')
fig.show()
