import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from scipy.signal import savgol_filter
import random
import multiprocessing
import sys
import math

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
    # n_cpu_max=multiprocessing.cpu_count() -2;
    n_cpu_max = 4
    print(n_cpu_max)
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


def gillespie_model(s0,t_obs_out,param):

    #--0--# Unpack parameters and species variables

    (p_LacI, p_TetR, p_CI, m_LacI, m_TetR, m_CI, m_GFP, p_GFP)=s0
    (k_m, k_m0, k_p, k_dm, k_dp, K, n, k_dGFP)=param

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

        types=["m_LacI_prod",
               "m_TetR_prod",
               "m_CI_prod",
               "m_LacI_loss",
               "m_TetR_loss",
               "m_CI_loss",
               "p_LacI_prod",
               "p_TetR_prod",
               "p_CI_prod",
               "p_LacI_loss",
               "p_TetR_loss",
               "p_CI_loss",
               "m_GFP_prod",
               "m_GFP_loss",
               "p_GFP_prod",
               "p_GFP_loss",
                ]

        #--1--#


        #--2--# Write rate expressions for each of the events

        rate_m_LacI_prod = k_m*K**n/ (K**n + p_CI**n)   + k_m0
        rate_m_TetR_prod = k_m*K**n/ (K**n + p_LacI**n) + k_m0
        rate_m_CI_prod   = k_m*K**n/ (K**n + p_TetR**n) + k_m0
        rate_m_GFP_prod  = k_m*K**n/ (K**n + p_TetR**n) + k_m0

        rate_p_LacI_prod = k_p*m_LacI
        rate_p_TetR_prod = k_p*m_TetR
        rate_p_CI_prod   = k_p*m_CI
        rate_p_GFP_prod  = k_p*m_GFP

        rate_m_LacI_loss = k_dm*m_LacI
        rate_m_TetR_loss = k_dm*m_TetR
        rate_m_CI_loss   = k_dm*m_CI
        rate_m_GFP_loss  = k_dm*m_GFP

        rate_p_LacI_loss = k_dp*p_LacI
        rate_p_TetR_loss = k_dp*p_TetR
        rate_p_CI_loss   = k_dp*p_CI
        rate_p_GFP_loss  = k_dGFP*p_GFP


        #--2--#



        #--3--# Store the rates into a list preserving the order of step 1.

        rates=[rate_m_LacI_prod,
               rate_m_TetR_prod,
               rate_m_CI_prod,
               rate_m_LacI_loss,
               rate_m_TetR_loss,
               rate_m_CI_loss,
               rate_p_LacI_prod,
               rate_p_TetR_prod,
               rate_p_CI_prod,
               rate_p_LacI_loss,
               rate_p_TetR_loss,
               rate_p_CI_loss,
               rate_m_GFP_prod,
               rate_m_GFP_loss,
               rate_p_GFP_prod,
               rate_p_GFP_loss,
            ]

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

        if event_type=="m_LacI_prod":
            m_LacI+=1;
            # pass
        elif event_type=="m_TetR_prod":
            m_TetR+=1;
            # pass
        elif event_type=="m_CI_prod":
            m_CI  += 1;
        elif event_type=="m_LacI_loss":
            m_LacI+= -1;

        elif event_type=="m_TetR_loss":
            m_TetR   += -1;

        elif event_type=="m_CI_loss":
            m_CI  += -1

        elif event_type=="p_LacI_prod":
            p_LacI  +=1;
        elif event_type=="p_TetR_prod":
            p_TetR  +=1;
        elif event_type=="p_CI_prod":
            p_CI    +=1;
        elif event_type=="p_LacI_loss":
            p_LacI  += -1;
        elif event_type=="p_TetR_loss":
            p_TetR  += -1;
        elif event_type=="p_CI_loss":
            p_CI    += -1
        elif event_type=="p_GFP_prod":
            p_GFP   += +1
        elif event_type=="p_GFP_loss":
            p_GFP   += -1
        elif event_type=="m_GFP_prod":
            m_GFP   +=  +1
        elif event_type=="m_GFP_loss":
            m_GFP   +=  -1;
        else:
            print "error unknown event type!!"

        #--4--#

        # store observation
        species=(p_LacI, p_TetR, p_CI, m_LacI, m_TetR, m_CI, m_GFP, p_GFP)

        t_obs.append(t)
        s_obs.append(species)

        # loops until time t exceeds t_final

    # loop has ended

    # before we return the results we must
    # resample the output to provide observations in accordance
    # with the t_obs passed to the function
    s_obs_out=resample_observations(t_obs,s_obs,t_obs_out)
    return np.array(s_obs_out)


# DEFINE INITIAL CONDITIONS AND PARAMETERS

#intitial condtions
p_LacI=2000
p_TetR=0
p_CI=0

m_LacI=0
m_TetR=0
m_CI=0

s=(p_LacI, p_TetR, p_CI, m_LacI, m_TetR, m_CI)

k_m=0.5
k_m0=5e-4
t_half_p=8.*60
t_half_m=2.*60
av_p_per_mRNA=20.0
K=40.0
n=2.1

k_dm=math.log(2.0)/t_half_m
k_dp=math.log(2.0)/t_half_p
k_dGFP=0.000128
t_av_p=1./k_dp
t_av_m=1./k_dm
k_p=av_p_per_mRNA/t_av_m


s0=(p_LacI, p_TetR, p_CI, m_LacI, m_TetR, m_CI, 0, 0)
param=    (k_m, k_m0, k_p, k_dm, k_dp, K, n, k_dGFP)


t_max=1000.*60.
# t_max=1000.*60.
t_obs=np.linspace(0,t_max,t_max/10+1) # store obs every 10 seconds
# t_obs = [0,1]
s_obs=gillespie_model(s0,t_obs,param)

# p_LacI, p_TetR, p_CI, m_LacI, m_TetR, m_CI, _ , _ = s_obs

p_LacI_obs = s_obs[:,0]
p_TetR_obs = s_obs[:,1]
p_CI_obs   = s_obs[:,2]
m_LacI_obs = s_obs[:,3]
m_TetR_obs = s_obs[:,4]
m_CI_obs   = s_obs[:,5]
m_GFP_obs   = s_obs[:,6]
p_GFP_obs   = s_obs[:,7]
