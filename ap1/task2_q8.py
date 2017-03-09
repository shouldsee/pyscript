import utils
from scipy.integrate import odeint
from utils import *
# %matplotlib inline

fss=[0.05,0.10,0.20,0.30,0.40];
fi=0.001;
siz=len(fss);
fig=plt.figure(figsize=[5,25]);

params=(5.,8.0);
N0=[(.199,0.001,0.800),('frac_S','frac_I','frac_R')];
s0=intobj(fNdot,N0,params);

for i in range(siz):
    fr=1-fss[i]-fi;
    s0.ss=(fss[i],fi,fr);
    s0.evo(np.linspace(0,100,1001))
    ax=plt.subplot(siz,1,i+1);
    
    ax.plot(s0.Ns[:,0],s0.Ns[:,1]);
        
#     s0.phase(ax,(0,1))
    ax.set_xlabel(s0.labels[0]);
    ax.set_ylabel(s0.labels[1]); 
    ax.set_title('Phase plot for with initial frac_s='+str(fss[i]))
    
    ax.set_ylim(0,0.15)
    ax.set_xlim(0,0.5)
fig.saveas('task2_q7.png')

fig=plt.figure(figsize=[5,5]);

ax=plt.subplot(1,1,1);
   
for i in range(siz):
    fr=1-fss[i]-fi;
    s0.ss=(fss[i],fi,fr);
    s0.evo(np.linspace(0,100,1001))
    ax.plot(fss[i],fi,'x');    
    ax.plot(s0.Ns[:,0],s0.Ns[:,1],label='initial_fs='+str(fss[i]));
        
#     s0.phase(ax,(0,1))
    ax.set_xlabel(s0.labels[0]);
    ax.set_ylabel(s0.labels[1]); 
    ax.set_title('Phase plots for different initial fs');
    ax.legend()
    
    ax.set_ylim(0,0.15)
    ax.set_xlim(0,0.5)
ax.plot([0.125,0.125],[0, 0.2])
# params=(5.,8.0)
# vfield(ax,fNdot,20,2/0.5,'len')
fig.savefig('task2_q8.png')