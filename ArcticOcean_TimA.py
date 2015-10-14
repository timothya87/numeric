
"""
Written by Tim Atkinson 
Mini Project 1

Purpose: to show the Arctic Ocean Temperatures profiles for a depth of 
200m

"""

import matplotlib.pyplot as plt
import numpy as np

plt.style.use("ggplot")

# constans
Td = 272.
Td1 = 271.
cp = 4.*10.**6.
beta = 0.8
alpha = 0.1
alb   = 0.1
h = 10
Amax = 1.*10.**-2.
Adep = 1.*10.**-4.
Adip = 1.5*10.**-3.


# create array of values for d
deltaD = 0.05
d = np.arange(0.+deltaD,200.,deltaD)
# and for I and DI/Dt
I = 100.*(1.-beta)*(1.-alb)*np.exp(-alpha*d)
dIdt = -alpha*100.*(1.-beta)*(1.-alb)*np.exp(-alpha*d)

# create all empty arrays for all variables
Ah  = np.empty_like(d)
dAh = np.empty_like(d)
T1  = np.empty_like(d)
T2  = np.empty_like(d)
T3  = np.empty_like(d)
b_array = np.empty_like(d)

# create all other arrays that vary for each depth,
# and differ the values for d<h, d>h.
for i in range(0,len(d)):
    
    if d[i] <= h:
        Ah[i] = Amax
        dAh[i] = 0. 
        T1[i] = 1.*(Ah[i]/(deltaD**2.))
        T2[i] = -2.*(Ah[i]/(deltaD**2.))
        T3[i] = 1.*(Ah[i]/(deltaD**2.))
        b_array[i] = (1./cp)*dIdt[i]

    else:
        Ah[i]   = Adep + ((Amax - Adep - Adip*(d[i]-h))*np.exp(-0.5*(d[i]-h)))
        dAh[i]  = -np.exp(-0.5*(d[i]-h))*(Adip + (0.5*(Amax-Adep-Adip*(d[i]-h))))
        T1[i]   = ((dAh[i]/(deltaD)) + (Ah[i]/deltaD**2.))
        T2[i]   = -(((dAh[i]/deltaD)) + (2.*Ah[i]/(deltaD**2.)))
        T3[i]   = Ah[i]/(deltaD**2.)
        b_array[i] = (1./cp)*dIdt[i]
    


def nmatrix(Ti1,Ti2,Ti3):
    """
    creates matrix of size n for problem 2 for 
    the varying depth steps for Ti+1, Ti, Ti-1
    """

    n = len(Ti1)
    a = np.zeros((n,n+2))
    r,c = np.indices(a.shape)
    
    a[r==c] = Ti1
    a[r==c-1] = Ti2
    a[r==c-2] = Ti3
    a = np.delete(a,0,1)
    a = np.delete(a,-1,1)
    a = np.array(a)
   
    return(a)
    


def smatrix(Ti1,Ti2,Ti3):
    """ 
    creates and array with the boundary conditions 
    substituted into the equations of size n for 
    depth steps for Ti+1, Ti, Ti-1
    """
    n = len(Ti1)
    a = np.zeros((n,n+2))
    r,c = np.indices(a.shape)
    a[r==c] = Ti1
    a[r==c-1] = Ti2
    a[r==c-2] = Ti3
    top = np.zeros(n+2)
    top[0] = 1.
    bot = np.zeros(n+2)
    bot[n+1] = 1.
    a = np.vstack([top,a,bot])
    
    return(a)


    

# append the boundary conditions to the b[] vector
b_array = np.append(Td,b_array)
b_array = np.append(b_array,Td1)

# create matrix using smatrix function
a = smatrix(T1,T2,T3)

# solve the matrix and save output as depths
temps = np.linalg.solve(a,b_array)
depths = np.append(0.,d)
depths = np.append(depths,200.)
    

# Plot it
fig,ax=plt.subplots(1,1,figsize=(5,5))
ax.plot(temps,depths)
ax.grid(True)
ax.set(xlabel="Temperature K",
       ylabel="Depth in meters",
       title="Arctic Ocean Temperatures")
plt.gca().invert_yaxis()


