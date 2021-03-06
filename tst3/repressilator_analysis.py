import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from scipy.signal import savgol_filter
import random
import multiprocessing
import sys
import math
import pickle

def save(object, filename, bin = 1):
	"""Saves an object to disk
	"""
	file = open(filename, 'wb')
	file.write(pickle.dumps(object, bin))
	file.close()


def load(filename):
	"""Loads an object from disk
	"""
	file = open(filename, 'rb')
	buffer = ""
	while 1:
		data = file.read()
		if data == "":
			break
		buffer += data
	object = pickle.loads(buffer)
	file.close()
	return object

# TEMPLATE TO RUN GILLESPIE SIMULATION OF MODEL
def peak_finder(timepoints, observations, threshold=50):
    maxima=[]
    minima=[]
    last_max=0
    last_min=0
    last_max_time=0
    last_min_time=0
    direction=1.0
    start_val=observations[0]
    for i in range(len(timepoints)):
        if start_val>start_val+threshold:
            direction=1.0
            break
        if start_val<start_val+threshold:
            direction=-1.0
            break
    for i in range(len(timepoints)):
        #print i, last_max, observations[i]
        if observations[i]*direction > last_max*direction:
            last_max=observations[i]
            last_max_time=timepoints[i]
        if observations[i]*direction < last_max*direction - threshold:
            if direction > 0:
                maxima.append((last_max_time,last_max))
            else:
                minima.append((last_max_time,last_max))
            last_max=observations[i]
            last_max_time=timepoints[i]
            direction*=-1.0
    return maxima

# CODE TO GENERATE THE DATA FILE

# s0=(p_LacI, p_TetR, p_CI, m_LacI, m_TetR, m_CI, 0, 0)
# param=    (k_m, k_m0, k_p, k_dm, k_dp, K, n, k_dGFP)

# t_transient=450.*60.
# s_obs_transient=gillespie_model(s0,[0,t_transient],param)
# s0=s_obs_transient[-1]

t_max=1000.*60.
t_obs=np.linspace(0,t_max,(t_max/10)+1)

# EXAMPLE CODE TO SIMULATE SINGLE RUN
# CODE TO SIMULATE MULTIPLE RUNS
# create a job pool to run the simulation multiple times
# add argument n_cpu=1 to disable multi core processing
#                           function name     num runs    arguments for function
# output=multiprocess(gillespie_model,   100,     fn_args=(s0,t_obs,param), n_cpu=4 )
# save(output,'repressilator_output.pickle')

output=load('repressilator_output.pickle')
print "# runs in output:",len(output)
print "# observations in each run",len(output[0])

# t_obs=np.linspace(0,t_max,t_max/10+1) # store obs every 10 seconds
# t_obs = [0,1]
# s_obs=gillespie_model(s0,t_obs,param)

# p_LacI, p_TetR, p_CI, m_LacI, m_TetR, m_CI, _ , _ = s_obs
plt.close("all")
fig1 = plt.figure(figsize=(8,6))
axL=fig1.add_subplot(1,1,1)
# axR=fig1.add_subplot(2,1,2)
axL.legend();
# axR.legend();

runnum=0;
for s_obs in output[:2]:
    runnum+=1;
    p_LacI_obs = s_obs[:,0]
    p_TetR_obs = s_obs[:,1]
    p_CI_obs   = s_obs[:,2]
    m_LacI_obs = s_obs[:,3]
    m_TetR_obs = s_obs[:,4]
    m_CI_obs   = s_obs[:,5]
    m_GFP_obs   = s_obs[:,6]
    p_GFP_obs   = s_obs[:,7]
    axL.plot(t_obs/60,p_LacI_obs,label = 'run'+str(runnum))

axL.set_xlabel('Time (min)')
axL.set_ylabel('LacI Proteins per cell')
axL.legend()
fig1.savefig('rep_runs.png')

plt.close("all")
fig2 = plt.figure(figsize=(8,6))
axL=fig2.add_subplot(1,1,1)
# axR=fig1.add_subplot(2,1,2)
axL.legend();

runnum=0;
for s_obs in output[:2]:
    runnum+=1;
    p_LacI_obs = s_obs[:,0]
    p_TetR_obs = s_obs[:,1]
    p_CI_obs   = s_obs[:,2]
    m_LacI_obs = s_obs[:,3]
    m_TetR_obs = s_obs[:,4]
    m_CI_obs   = s_obs[:,5]
    m_GFP_obs   = s_obs[:,6]
    p_GFP_obs   = s_obs[:,7]
    axL.plot(t_obs/60,p_LacI_obs,label = 'run'+str(runnum))

avg = output.mean(axis = 0);
s_obs = avg;
p_LacI_obs = s_obs [:,0]

axL.plot(t_obs/60,p_LacI_obs,label = 'avg_run')
axL.legend()
fig2.savefig('rep_runs(updated).png')


