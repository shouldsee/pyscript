def fNdot(N,t,params):
#     tr=T, r0=r0;
    (s,i,r)=N;
    T,r0=params;
    ds=0.      -r0*s*i/T;
    di=r0*s*i/T-i/T;
    dr=i/T;
    return (ds,di,dr)
params=(5.,8.0);
a=.2;
b=.001
N0=[(.199,0.001,0.800),('frac_S','frac_I','frac_R')];
s1=intobj(fNdot,N0,params);
s1.evo(np.linspace(0,100,1001));
fig2=plt.figure(figsize=[10,9])

ax=plt.subplot(2,2,2)
for i in range(3):
    s1.line(ax,i);
ax.legend()
ax.set_xlabel('time(days)')
ax.set_ylabel('f_s')
ax.set_title('all population vs time for SIR model')
ax.legend();

ax=plt.subplot(2,2,1)
for i in range(1):
    s1.line(ax,i);
ax.legend()
ax.set_xlabel('time(days)')
ax.set_ylabel('f_s')
ax.set_title('susceptible population vs time for SIR model')
ax.legend();

ax2=plt.subplot(2,2,3)
s1.line(ax2,1);
ax2.set_xlabel('time(days)')
ax2.set_ylabel('f_i')
ax2.set_title('Infected population vs time for SIR model')
ax2.legend();
ax2.legend(loc=1)
fig2.savefig('task2_q4.png')