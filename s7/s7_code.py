from myreport import html_report
import random
from matplotlib import pyplot as plt
import numpy as np
plt.style.use('ggplot')

# PUT ALL FUNCTIONS HERE

def draw_ball(G,R):
    pG=float(G)/float(G+R)
    rand_num=random.random()
    if rand_num<pG:
        return 'G'
    else:
        return 'R'

def simulate_system():
    #balls
    G=10.
    R=5.
    T=G+R

    #time
    t=0

    G_record=[]
    R_record=[]
    t_record=[]

    G_record.append(G)
    R_record.append(R)
    t_record.append(t)

    while t<=1000:
   	# increment time
   	t=t+10.#/100.
        if G+R!=0:
       	    # select ball from bag
       	    ball=draw_ball(G,R)
       	    # update number of balls
            if ball=='G':
                G=G-1
                R=R+1
   	# record system state
   	G_record.append(G)
   	R_record.append(R)
   	t_record.append(t)
    return (G_record, R_record, t_record)

#results= simulate_system()
#G_record, R_record, t_record = results
#plt.plot(t_record, G_record, 'g-')
#plt.plot(t_record, R_record, 'r-')
#plt.show()

my_runs=[]
t_record=None
for i in range(200):
	results= simulate_system()
	G_record, R_record, t_record = results
	my_runs.append(G_record)

print "length t_record:",len(t_record)
print "length my_runs:",len(my_runs)
print "length my_runs[0]:",len(my_runs[0])

plt.close('all')

# convert to array and transpose
my_runs=np.array(my_runs)
my_runs=np.transpose(my_runs)

fig1=plt.figure()
ax=fig1.add_subplot(1,1,1)
ax.plot(t_record, my_runs, '-')

av_run=np.average(my_runs,axis=1)
# axis=0 averages colums
# axis=1 averages rows
ax.plot(t_record, av_run, 'k-', lw=2)
fig1.show()

fig2=plt.figure()
ax=fig2.add_subplot(1,1,1)
ax.hist(my_runs[10,:], bins=np.arange(-0.5,11.5,1.0))
fig2.show()

# now write results
myreport=html_report("s7_starting_code.html")
myreport.add_subheading('Decay trajectories')
myreport.add_figure(fig1)
myreport.add_subheading('Distribution of G at t=100 (200 runs)')
myreport.add_figure(fig2)

myreport.add_subheading('Python Code')
myreport.add_source(__file__)
myreport.write()
myreport.view()
