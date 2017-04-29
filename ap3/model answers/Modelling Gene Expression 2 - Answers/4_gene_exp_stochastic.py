import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
import random
import multiprocessing
import sys

# MODEL FOR STOCHASTIC GENE EXPRESSION

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


###################
# Helper functions (Do not change!)

def resample_observations(t_in_obs,s_in_obs,t_out_obs):
    s_out_obs=[]
    i=0
    ti=t_in_obs[i]
    si=s_in_obs[i]
    j=0
    while j<len(t_out_obs):
        T=t_out_obs[j]
        while ti<T and i<len(t_in_obs)-1:
            si=s_in_obs[i]
            i+=1
            ti=t_in_obs[i]
        s_out_obs.append(si)
        j+=1
    return s_out_obs


def gen_next_event_time(rate):
    t=random.expovariate(rate)
    return t


def random_choice_from_pdf(pdf):
    cdf=[]
    cumulative_p=0
    for p in pdf:
        cumulative_p+=p
        cdf.append(cumulative_p)
    rand=random.random()

    for i in range(len(cdf)):
        if rand<cdf[i]:
            return i
    # last cdf should be 1.0 so the following should never happen!
    print "Error generating choice, check PDF"
    return None

def multiprocess(fn, n_runs, fn_args, n_cpu=None):
    n_cpu_max=multiprocessing.cpu_count()
    if n_cpu is None or n_cpu > n_cpu_max:
        n_cpu=n_cpu_max
        print "Using",n_cpu,"processing thread(s)"
    queue=[]
    for i in range(n_runs):
        queue.append({'id':i,'fn':fn,'args':fn_args,'total_jobs':n_runs})
    if n_cpu!=1:
        p=multiprocessing.Pool(n_cpu)
        result=p.map(run_multiproc,queue, chunksize=1)
        p.close()
        p.join()
    else:
        result=[]
        for job in queue:
            result.append(run_multiproc(job))
    return result

def run_multiproc(job_dict):
    fn=job_dict['fn']
    s =  'Progress: job {:3d}/{:3d} running'.format(job_dict['id']+1,job_dict['total_jobs'])
    sys.stdout.write('\b'*len(s))
    sys.stdout.write(s)
    sys.stdout.flush()
    result=fn(*job_dict['args'])
    return result

###################


def gillespie_gene_expression(s0,t_obs_out,param):

    #--0--# Unpack parameters and species variables

    k_on, k_off, k_t, k_p, kd_mRNA, kd_p, K, n = param
    g_on,g_off,mRNA,p = s0

    #--0--#

    # create arrays for output
    s_obs=[]
    t_obs=[]

    # read in start time and end time
    t_init=t_obs_out[0]
    t_final=t_obs_out[-1]

    t=t_init
    t_obs.append(t)
    s_obs.append(s0)

    while t < t_final:

        #--1--# Write labels for each event type here.

        types=['on','off','t','p','dmRNA','dp']

        #--1--#



        #--2--# Write rate expressions for each of the events

        rate_on = g_off*k_on
        rate_off = g_on*k_off
        rate_t = g_on*k_t*K**n/(K**n+p**n)
        rate_p = mRNA*k_p
        rate_dmRNA = mRNA*kd_mRNA
        rate_dp = p*kd_p

        #--2--#



        #--3--# Store the rates into a list preserving the order of step 1.

        rates=[rate_on,rate_off,rate_t,rate_p,rate_dmRNA,rate_dp]

        #--3--#


        #-- Do not edit below --#

        ## CARRY OUT GILLESPIE ALGORITHM TO STEP FORWARD TO NEXT EVENT
        ## AND UPDATE SYSTEM STATE ACCORDING TO EVENT TYPE

        # calc total reaction rate
        rate_all_events=sum(rates)

        # if rate of events is zero break from loop
        # e.g. when all reactants used up
        if rate_all_events==0:
            break

        # generate the time until the next event
        # in accordance with rate_all_events
        next_event=gen_next_event_time(rate_all_events)

        # calc PDF for event type
        # in accordance with relative rates
        pdf=[]
        for event_rate in rates:
            p_event = event_rate/sum(rates)
            pdf.append(p_event)

        rand_i =  random_choice_from_pdf(pdf)
        event_type=types[rand_i]

        # increment time and number of molecules
        # according to event type
        t=t+next_event

        #-----------------------------------#



        ## ALGORITHM HAS INCREMENTED TIME AND SELECTED NEXT EVENT
        ## WE NOW NEED TO UPDATE OUR SYSTEM ACCORDING TO THE EVENT
        ## TYPE STORED IN VARIABLE event_type



        #--4--# Complete the if-elif-else commands to update the system
              # according to event type

        if event_type=="on":
            g_on+=1
            g_off+=-1
        elif event_type=="off":
            g_off+=1
            g_on+=-1
        elif event_type=="t":
            mRNA+=1
        elif event_type=="p":
            p+=1
        elif event_type=="dmRNA":
            mRNA+=-1
        elif event_type=="dp":
            p+=-1
        else:
            print "error unknown event type!!"

        #--4--#

        # store observation
        species=(g_on,g_off,mRNA,p)

        t_obs.append(t)
        s_obs.append(species)

        # loops until time t exceeds t_final

    # loop has ended

    # before we return the results we must
    # resample the output to provide observations in accordance
    # with the t_obs passed to the function
    s_obs_out=resample_observations(t_obs,s_obs,t_obs_out)
    return np.array(s_obs_out)


#intitial condtions
g_off=0.0
g_on=1.0
mRNA=0
p=0


k_on=1./(10.)
k_off=1./(2.)
k_t=1./20.
k_p=1./2.0
kd_mRNA=1./(8.*60.)
kd_p=1./(50.*60.)
K=4000.
n=2

s0=(g_on,g_off,mRNA,p)
param=(k_on, k_off, k_t, k_p, kd_mRNA, kd_p,K,n)

t_max=30000
t_obs=np.linspace(0,t_max,t_max+1)

mRNA_obs=[]
#t_final=30000
#t_obs=np.linspace(0,t_final,t_final+1.0)


# EXAMPLE CODE TO SIMULATE SINGLE RUN
#s_obs=gillespie_gene_expression(s0,t_obs,param)

# CODE TO SIMULATE MULTIPLE RUNS
# create a job pool to run the simulation multiple times
#
# add argument n_cpu=1 to disable multi core processing
# if you experience crashes
#
#                           function name     num runs    arguments for function
output=multiprocess(gillespie_gene_expression, 500, fn_args=(s0,t_obs,param))

g_on_obs=[]
mRNA_obs=[]
p_obs=[]

for run in output:
    g_on_obs.append(run[:,0])
    mRNA_obs.append(run[:,2])
    p_obs.append(run[:,3])

# convert to np arrays to allow us to manipulate them
g_on_obs=np.array(g_on_obs)
mRNA_obs=np.array(mRNA_obs)
p_obs=np.array(p_obs)

# transpose so it is arranged with:
# observations points as rows (first index)
# run number as columns (second index)
g_on_obs=np.transpose(g_on_obs)
mRNA_obs=np.transpose(mRNA_obs)
p_obs=np.transpose(p_obs)

plt.close("all")
fig0=plt.figure()

# add "super title" over plot set
fig0.suptitle('Stochastic Gene Expression Model')
ax_G=fig0.add_subplot(3,1,1)
ax_mRNA=fig0.add_subplot(3,1,2)
ax_P=fig0.add_subplot(3,1,3)

ax_G.set_ylabel('fraction $g_{on}$')
ax_mRNA.set_ylabel('#mRNA copies')
ax_P.set_ylabel('#Protein copies')
ax_P.set_xlabel('time (s)')

# use array slice [:,0]
# which plots all time bservations for run 0
ax_G.plot(t_obs, g_on_obs[:,0], 'b-')
ax_mRNA.plot(t_obs, mRNA_obs[:,0], 'r-')
ax_P.plot(t_obs, p_obs[:,0], 'g-')

# adjusts matplotloib to autoscale margins with
# 0.1 i.e 10% head room at the top of the y-axis
ax_G.margins(y=0.1)
ax_G.set_ylim(bottom=0)
ax_P.margins(y=0.1)
ax_P.set_ylim(bottom=0)
ax_mRNA.margins(y=0.1)
ax_mRNA.set_ylim(bottom=0)

fig0.show()


# we want to select all run vals from the last observation point (time=3000)
# observations are first index, runs second index
p_end=p_obs[-1,:]
mRNA_end=mRNA_obs[-1,:]


fig1=plt.figure()
ax_P=fig1.add_subplot(1,1,1)
ax_P.hist(p_end, bins=np.arange(0,15.1,1))
ax_P.set_title("Histogram of P copy number at t=30000, (500 runs)")
ax_P.set_xlabel('#  proteins')
fig1.show()

fig2=plt.figure()
ax_mRNA=fig2.add_subplot(1,1,1)
ax_mRNA.hist(mRNA_end, bins=10)
ax_mRNA.set_title("Histogram of mRNA copy number at t=30000, (500 runs)")
ax_mRNA.set_xlabel('#  mRNA copies')
fig2.show()
mean_mRNA=np.mean(mRNA_end)
mean_P=np.mean(p_end)

stdev_mRNA=np.std(mRNA_end)
stdev_P=np.std(p_end)

print "Mean mRNA copies: {:.3g}, stdev:{:.3g}".format(mean_mRNA,stdev_mRNA)
print "Mean Protein copies: {:.4g}, stdev:{:.4g}".format(mean_P,stdev_P)

CV_mRNA=stdev_mRNA/mean_mRNA
CV_P=stdev_P/mean_P
print "CV mRNA: {:.3g}".format(CV_mRNA)
print "CV P: {:.3g}".format(CV_P)
