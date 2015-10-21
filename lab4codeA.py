# Written by Tim Atkinson 10526135
# Atsc 409

# import packages 
import json
import numpy as np
from matplotlib import pyplot as plt
from collections import namedtuple


# wrote out functions in script due to import errors in previous attempts
def read_init(filename):
  with open(filename,'r') as f:
      init_dict=json.load(f)
  #
  #convert dictionary to namedtuple
  #
  initvals=namedtuple('initvals','dt c1 c2 t_beg t_end yinitial comment plot_title')
  theCoeff=initvals(**init_dict)
  return theCoeff

# Heun's method
#*****************************************
def heun(coeff, y, derivs):
  ynew = y + coeff.dt*derivs(coeff,y + ((2/3) * coeff.dt * derivs(coeff,y)))
  return ynew
#*****************************************
  
# midpoint method
def midpoint(coeff, y, derivs):
  ynew = y + coeff.dt*derivs(coeff,y + ((1/2) * coeff.dt * derivs(coeff,y)))
  return ynew  

# dervis function to solve differential equation
def derivs(coeff, y):
  f=np.empty_like(y) #create a 2 element vector to hold the derivitive
  f[0]=y[1]
  f[1]= -1.*coeff.c1*y[1] - coeff.c2*y[0]
  return f
#
# first make sure we have an input file in this directory
#
initialVals={'yinitial': [0.,1.],'t_beg':0.,'t_end':40.,'dt':0.1,'c1':(1/4),'c2':(3/4)}
initialVals['comment'] = 'written Sep. 29, 2015'
initialVals['plot_title'] = 'Damped Oscilator'

infile='run_1.json'
with open(infile,'w') as f:
      f.write(json.dumps(initialVals,indent=4))
#
#  now read the initial information into a namedtuple coeff
#
        
infile='run_1.json'
coeff=read_init(infile)

#
# integrate and save the result in savedata
#
# Plot the damped Oscilator 
time=np.arange(coeff.t_beg,coeff.t_end,coeff.dt)
y=coeff.yinitial
nsteps=len(time) 
savedata=np.empty([nsteps],np.float64)
for i in range(nsteps):
    y=heun(coeff,y,derivs)
    savedata[i]=y[0]

theFig,theAx=plt.subplots(1,1,figsize=(8,5))
theAx.plot(time,savedata,'o-')
theAx.set_title(coeff.plot_title)
theAx.set_xlabel('time (seconds)')
theAx.set_ylabel('y0')       
