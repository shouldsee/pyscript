from reportlib import html_report
import numpy as np
from matplotlib import pyplot as plt
from scipy.integrate import odeint
from matplotlib.widgets import Slider, Button, RadioButtons

## OUR ODE MODEL FUNCTION

def sdot(s,t,p):
    H,L=s
    k1,k2,k3,eff=p

    dH = k1*H - k2*H*L
    dL = eff*k2*H*L - k3*L

    dS=(dH,dL)
    return dS

H0=100
L0=1
k1=0.7
k2=0.3
k3=0.5
eff=0.02

## RUN MODEL
s0=(H0,L0)
p=(k1,k2,k3,eff)
t_max=20
t_obs=np.linspace(0,t_max,1001)
s_obs=odeint(sdot,s0,t_obs,args=(p,))

LV_report=html_report("LV.html")

## PUT SOME TEXT INTO OUR REPORT
LV_report.add_heading("The Lotka-Volterra Model")
LV_report.add_subheading("Equations:")
LV_report.add_code("dH/dt = k1 * H        -  k2 * H * L")
LV_report.add_code("dL/dt = eff * k2 . H  -  k3 * L")

## CREATE A FIGURE FOR OUR REPORT

## BECAUSE FIGURE IS GOING IN REPORT
## WE CREATE IT LIKE THIS:
fig = LV_report.init_figure()
ax=fig.add_subplot(2,1,1)
ax2=fig.add_subplot(2,1,2)
ax.set_title('Hare')
ax2.set_title('Lynx')
ax.plot(t_obs, s_obs[:,0], 'b-')
ax2.plot(t_obs, s_obs[:,1], 'g-')

## INSERT OUR FIGURE INTO THE REPORT
LV_report.add_subheading('Output of Model')
LV_report.add_text('Parameters: k1,k2,k3,eff=0.7,0.3,0.5,0.02')
LV_report.add_text('Initial conditions: H0,L0=100,1')
LV_report.add_figure(fig)

# INSERT SOURCE CODE INTO REPORT
LV_report.add_subheading('Python Code')
LV_report.add_source(__file__)

# WRITE REPORT TO FILE
LV_report.write()

# OPEN FILE IN BROWSER
LV_report.view()
