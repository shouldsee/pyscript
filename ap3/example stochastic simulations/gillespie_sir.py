import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
import random
import multiprocessing
import sys

# TEMPLATE TO RUN GILLESPIE SIMULATION OF MODEL

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
    if n_cpu!=1 and len(queue)>1:
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


def gillespie_sir(s0,t_obs_out,param):

    #--0--# Unpack parameters and species variables

    ki, kr = param
    S,I,R = s0

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

        types=['infection','recovery']

        #--1--#


        #--2--# Write rate expressions for each of the events

        rate_infection = ki*S*I
        rate_recovery = kr*I

        #--2--#



        #--3--# Store the rates into a list preserving the order of step 1.

        rates=[rate_infection, rate_recovery]

        #--3--#


        #-- Do not edit the section below --#

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

        if event_type=='infection':
            S+=-1
            I+=1
        elif event_type=='recovery':
            I+=-1
            R+=1
        else:
            print "error unknown event type!!"

        #--4--#

        # store observation
        species=(S,I,R)

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
S0=1990
I0=10
R0=8000
N=S0+I0+R0

#parameters
ki=10.0/5.0*1./N
kr=1./5.

s0=(S0,I0,R0)
param=(ki,kr)

t_max=150
t_obs=np.linspace(0,t_max,t_max+1)

# EXAMPLE CODE TO SIMULATE SINGLE RUN
#s_obs=gillespie_gene_expression(s0,t_obs,param)

# CODE TO SIMULATE MULTIPLE RUNS
# create a job pool to run the simulation multiple times
# add argument n_cpu=1 to disable multi core processing
#                   function na  num runs    arguments for function
output=multiprocess(gillespie_sir,   1,     fn_args=(s0,t_obs,param) )

# example code for simple plot of observations

s_obs=output[0]

plt.close("all")

fig=plt.figure()

ax=fig.add_subplot(1,1,1)

ax.plot(t_obs, s_obs[:,1])

fig.show()
