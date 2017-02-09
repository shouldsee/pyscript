import math
import numpy as np
from matplotlib import pyplot as plt
from functools import reduce
def get_factorial(n):
	if n<0 or int(n)!=n:
		print("factorial input invalid");
		
	facts=range(1,n+1);
	if not facts:
		nfact=1;
	else:
		nfact=reduce( lambda x,y:x*y, facts);
	return nfact;


def taylor_expansion_exp(dt,n):
	term=dt**n/get_factorial(n)
	return term

# now write function that sums a set of terms (i=0 to i=n)

def sum_taylor_expansion_exp(dt, max_n):
	summed_terms=0
	for ni in range(0, max_n+1):
		summed_terms+=taylor_expansion_exp(dt,ni)
	return summed_terms

# now write master function to plot Taylor approximation

# MASTER FUNCTION

def plot_taylor_approximation_exp(n, tmax):

	# plot the taylor expansion of the exponential function up to term n
	#
	# showComponents = True or False
	# turn on/off additional plot of each term in the series
	#
	# showError = True / False
	# turn on/off additional plot showing relative error

	tmin=0
	tmax=tmax

	dt=np.linspace(tmin,tmax,400)

	t=[]
	y_taylor=[]
	y_exact=[]

	for dti in dt:
		yi_taylor = sum_taylor_expansion_exp(dti,n)
		yi_exact=math.exp(dti)

		t.append(dti)
		y_taylor.append(yi_taylor)
		y_exact.append(yi_exact)

	plt.close("all")
	fig=plt.figure(figsize=[6,10])
	ax=fig.add_subplot(3,1,1)
	ax.set_title("Plot of Taylor approximation (Exact function shown as dashed line)",fontsize='small')
	ax.plot(t,y_exact,'-k',ls='dashed',label="y=e^t")
	ax.plot(t,y_taylor,'-r',label="taylor approx to term {}".format(n))
	ax.set_xlim([tmin,tmax])
	ax.set_ylim([0,1.1*math.exp(tmax)])
	ax.legend(loc='lower right',fontsize='small')

	ax2=fig.add_subplot(3,1,2)
	ax2.set_title("Components of Taylor series",fontsize='small')
	components={}
	for i in range(n+1):
		components[i]=[]
	for dti in dt:
		for i in range(n+1):
			components[i].append(   taylor_expansion_exp(dti,i) )
		for i in range(n+1):
			ax2.plot(t,components[i],ls='dotted',label="{}th term".format(i))
	ax2.set_ylim([0,1.2])
	ax2.legend(loc='upper left',fontsize='small')
	ax3=fig.add_subplot(3,1,3)
	ax3.set_title("Plot of relative error (line is 1% level)",fontsize='small')
	rel_error=[]
	for i in range(0,len(t)):
		rel_error.append((y_exact[i]-y_taylor[i])/y_exact[i])
		ax3.plot(t,rel_error,label="rel err".format(n))
		ax3.set_ylim([0,0.1])
		ax3.plot([0,tmax],[0.01,0.01],ls='dotted',label="1% error")
		ax3.legend(loc='upper left',fontsize='small')

	fig.show()
	return True

# code to call master function
plot_taylor_approximation_exp(n=3, tmax=1.5)
