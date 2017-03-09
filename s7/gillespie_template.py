from myreport import html_report
import random
from matplotlib import pyplot as plt
import numpy as np
plt.style.use('ggplot')

def random_choice_from_pdf(pdf):
    cdf=[]
    cumulative_p=0
    for p in pdf:
        cumulative_p+=p
        cdf.append(cumulative_p)
    rand=random.random()
    for i in range(len(cdf)):
        if rand<=cdf[i]:
            return i
    # last cdf should be 1.0 so the following should never happen!
    print "Error generating choice, check PDF"
    return None


def simulate_system():
    #balls
    G=10.
    R=5.
    N=G+R
    k=1.0/10.*1.0/N

    #time
    t=0

    G_record=[]
    R_record=[]
    t_record=[]

    G_record.append(G)
    R_record.append(R)
    t_record.append(t)

    while t<=1000:

        # event types:
        # GE and RE

        # calculate rates
        rateGE=...
        rateRE=...

        # calculate total rate
        rateT=...

        P_GE=...
        P_RE=...

   	# generate time interval
   	interval=random.expovariate(rateT)

   	# increment time
   	t=t+...

   	# generate PDF
   	myPDF = ( ..., ... )

        # select event type
       	index=random_choice_from_pdf(myPDF )
        # index=0 corresponds to item 0 in PDF: GE event
        # index=1 corresponds to item 1 in PDF: RE event

        event=None

        if index==0:
            event=...
        elif index==1:
            event=...

       	# update number of balls
        if event=='GE':
            # update G and R numbers after GE event
            ...
        elif event=='RE':
            # update G and R numbers after RE event
            ...
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

G_runs=[]
t_runs=[]
t_record=None
half_life_measurements=[]
for i in range(1):
	results= simulate_system()
	G_record, R_record, t_record = results
	G_runs.append(G_record)
	t_runs.append(t_record)

print "length G_runs:",len(G_runs)
print "length t_runs:",len(t_runs)

print "length t_runs[0]:",len(t_runs[0])
print "length G_runs[0]:",len(G_runs[0])

plt.close('all')

fig1=plt.figure()
ax=fig1.add_subplot(1,1,1)
n_runs=len(t_runs)
for i in range(n_runs):
    ax.plot(t_runs[i], G_runs[i], '-')
fig1.show()

# now write results
myreport=html_report("s7_Gillespie.html")
myreport.add_subheading('Decay trajectories')
myreport.add_figure(fig1)

myreport.add_subheading('Python Code')
myreport.add_source(__file__)
myreport.write()
myreport.view()
