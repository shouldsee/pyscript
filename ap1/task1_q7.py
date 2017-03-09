
# coding: utf-8

# In[1]:

from task1_q5 import *


# In[25]:

s1.ss=(1990,10,8000);
s1.evo(np.linspace(0,100,1001));
fig2=plt.figure(figsize=[5,9])
ax=plt.subplot(2,1,1)
for i in range(3):
    s1.line(ax,i);
ax.legend()
ax.set_xlabel('time(days)')
ax.set_ylabel('N')
ax.set_title('All populations vs time for SIR model')
ax.legend();

ax2=plt.subplot(2,1,2)
s1.line(ax2,1);
ax2.set_xlabel('time(days)')
ax2.set_ylabel('I')
ax2.set_title('Infected vs time for SIR model')
ax2.legend();
ax2.legend(loc=1)
fig2.savefig('task1_q7.png')
print('total infection is ',s1.ss[0]-s1.Ns[-1,0])
    
    


# In[ ]:



