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
