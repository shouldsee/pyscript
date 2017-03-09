from myreport import html_report
import random
from matplotlib import pyplot as plt

# PUT ALL FUNCTIONS HERE
#def draw_ball():

#def draw_ball2(G,R):

def simulate_system():
    #balls
    G=10
    R=5
    T=G+R

    #time
    t=0

    G_record=[]
    R_record=[]
    t_record=[]

    G_record.append(G)
    R_record.append(R)
    t_record.append(t)

    while t<1000:
   	# increment time
   	t=t+10

       	# select ball from bag

       	# update number of balls

   	# record system state
   	G_record.append(G)
   	R_record.append(R)
   	t_record.append(t)
    return (G_record, R_record, t_record)

results= simulate_system()
G_record, R_record, t_record = results
plt.plot(t_record, G_record, 'g-')
plt.plot(t_record, R_record, 'r-')
plt.show()
