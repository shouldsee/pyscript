def Ndot(N,t,params):
    N1,N2=N;
    k1,k2,k3,eff=params;
    N1dot=-k3*N1+eff*k2*N1*N2;
    N2dot=k1*N2-k2*N1*N2;
    return (N1dot,N2dot);

N0=(25,1250);
params=(0.87,0.08,0.33,0.0025);
k1,k2,k3,eff=params;
# numstep
t_max=5000.;
t_stepsize=0.5;
t_num=round(t_max/t_stepsize)
ts=np.linspace(0,t_max,t_num);

def nullclines(ax):
    ax.plot([N1f,N1f],list(ax.get_ylim()),'r--',lw=1,label='null_1')
    ax.plot([0,0],list(ax.get_ylim()),'b--',lw=1,label='null_2')
    ax.plot(list(ax.get_xlim()),[N2f,N2f],'b--',lw=1)
    ax.plot(list(ax.get_xlim()),[0,0],'r--',lw=1)
    
# Use np.arange(t_min,t_max,t_step) instead
# def ls(t_min,t_max,t_step):
#     return np.linspace(t_min,t_max,int((t_max-t_min)/t_step))

def vfield(axis,Ndot=Ndot,num=20,scale=1/0.0015,norm_method='none'):
    xlim=list(axis.get_xlim());
    ylim=list(axis.get_ylim());
    x=np.linspace(xlim[0],xlim[1],num);
    y=np.linspace(ylim[0],ylim[1],num);
    xs,ys=np.meshgrid(x,y);
    us,vs=Ndot((xs,ys),0,params);
    
    ls=(us**2+vs**2)**0.5;        
    if norm_method=='log':
        lls=np.min(np.log(ls),0);
        us=us/ls*lls;
        vs=vs/ls*lls;
    elif norm_method=='len':
        ls=(us**2+vs**2)**0.5;
        us=us/ls;
        vs=vs/ls;
    axis.quiver(xs,ys,us,vs,color='b',scale_units='inches',angles='xy',pivot='mid',scale=scale)

plt.close("all");
fig=plt.figure(figsize=[20,10])
ax=plt.subplot(1,2,1)
ax2=plt.subplot(1,2,2)
Ns=odeint(Ndot,N0,np.arange(0,10,0.5),args=(params,))
x=Ns[:,0];
y=Ns[:,1];
ax.plot(x,y,lw=0.5,markersize=0.5)
ax.quiver(x[:-1], y[:-1], x[1:]-x[:-1], y[1:]-y[:-1],angles='xy',scale_units='xy',scale=1);
vfield(ax,Ndot,20,5/0.5,'log')
# nullclines(ax)

ax2.plot(x,y,lw=0.5,markersize=0.5)
ax2.quiver(x[:-1], y[:-1], x[1:]-x[:-1], y[1:]-y[:-1],angles='xy',scale_units='xy',scale=1);
vfield(ax2,Ndot,20,2/0.5,'len')
# nullclines(ax2)


