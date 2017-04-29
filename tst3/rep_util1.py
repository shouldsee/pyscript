import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
import math,copy,random
import sys
#%matplotlib inline


global IPTG

cmd_init_name='''
varname=['p_LacI','p_TetR','p_CI','m_LacI','m_TetR','m_CI','m_GFP','p_GFP'];
cstname=['k_m','k_m0','k_p','k_dm','k_dp','K','n','k_dGFP','IPTG','K_IPTG','K_TetR','cpnum'];
'''
cmd_name2dict='''
vardict={x:varname.index(x) for x in varname};
cstdict={x:cstname.index(x) for x in cstname};
'''
make_init_cmd = lambda varname,cstname:'s0=(%s);param = (%s);'% (','.join(varname),','.join(cstname));
make_unpack_cmd = lambda varname,cstname:'(%s)=s;(%s)=param;'% (','.join(varname),','.join(cstname));
make_pack_cmd = lambda varname,cstname: 'sdot=(%s);'%(','.join(['d'+x for x in varname]))
make_unpackobs_cmd = lambda vardict:''.join( '%s_obs = s_obs[:,%s];' % (key,val) for key,val in vardict.items());

cmd_default_incond='''
p_LacI=0
p_TetR=0
p_CI=0
p_GFP=0

m_LacI=5
m_TetR=0
m_CI=0
m_GFP=0
'''

cmd_default_param='''

cpnum0=5.;
k_m=0.5/cpnum0;
k_m0=5e-4/cpnum0;

K=40.0
K_TetR=40.0;
n=2.1

k_dm=0.00577622650467
k_dp=0.00144405662617
k_dGFP=0.000128360588993
k_p=0.115524530093
IPTG=0.;
K_IPTG=6.E6;
cpnum = 5;
'''
exec(cmd_default_incond)
exec(cmd_default_param)

exec(cmd_init_name)
exec(cmd_name2dict)


cmd_unpack_all=make_unpack_cmd(varname,cstname);
cmd_pack_all=make_pack_cmd(varname,cstname);
cmd_unpack_obs=make_unpackobs_cmd(vardict)
cmd_unpack_obs
# cmd_pack_all




def fft_period(ys,threshold=20, debug=1):
    X = np.array(ys)
    N=len(X)
    W    = np.fft.fft(X)/X.size
    W[0]=0;
    freq = np.fft.fftfreq(N,1)

#     threshold = 20
    absW=abs(W)
    try:
        dur = np.where(absW>threshold)[0];
        Wdur=absW[dur];
        idx = dur[np.where(Wdur == max(Wdur))[0]];
    #     idx = np.where(abs(W)>threshold)[0][-1]
#         axRa.plot(Wdur)
        max_f = abs(freq[idx])
    except:
        if debug:
            print('errored or no significant period');
        max_f = np.array(0);

    period = (1/max_f).flat[0];
    if debug:
        fig=plt.figure()
        ax = fig.add_subplot(1,1,1)
        ax.plot(1/freq[:N/2],absW[:N/2])
        print "Period estimate: ", period
    return period


def sdot(s,t,param):

#     (p_LacI, p_TetR, p_CI, m_LacI, m_TetR, m_CI, m_GFP, p_GFP)=s
#     (k_m, k_m0, k_p, k_dm, k_dp, K, n, k_dGFP, IPTG,fix_conc)=param
    exec(cmd_unpack_all)

#     if not fix_list==None:
#     print(fix_list)
    
#     for var_ind in fix_list:
#         cmd=varname[var_ind] + ' = (fix_conc[var_ind])';
#         exec(cmd);
#         print(cmd)
#         print(var_ind)
        
    rate_CIO = (k_m / (1. + (p_CI/K)**n) +
                k_m0);
    rate_LacO = (k_m*K**n / (K**n + p_LacI**n) +
                 k_m0);
    rate_LacO = (k_m/ (1. + (p_LacI/K * ( 1. /(1.+IPTG/K_IPTG))**2 )**n) +
                 k_m0);

    rate_TetO = (k_m / (1 + (p_TetR/K_TetR)**n) +                      
                 k_m0);
    
    rate_m_TetR_prod =  rate_LacO*cpnum;
    rate_m_LacI_prod = rate_CIO*cpnum;
    rate_m_CI_prod   = rate_TetO*cpnum;
    rate_m_GFP_prod  = 0

#     rate_m_TetR_prod = k_m*K**n / (K**n + p_LacI**n) + k_m0
#     rate_m_LacI_prod = k_m*K**n / (K**n + p_CI**n)   + k_m0
#     rate_m_CI_prod   = k_m*K**n / (K**n + p_TetR**n) + k_m0
    
    rate_p_LacI_prod = k_p*(m_LacI)
    rate_p_TetR_prod = k_p*(m_TetR)
    rate_p_CI_prod   = k_p*(m_CI)
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

    exec(cmd_pack_all)
    return sdot


def analyse(ys):
    period = fft_period(ys,debug=0)/nt;
    amp=ys.max()-ys.min()
    avg=ys[-int(period*nt):].mean();
    print('log_par = {}, period estimate: {:f}, amplitude:{:f}, mean_LacI:{:f}'.format(str(log_par),period,amp,avg))
    return((period,amp,avg))