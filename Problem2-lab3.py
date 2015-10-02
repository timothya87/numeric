# -*- coding: utf-8 -*-
"""
Written by Tim Atkinson 
Problem 2 a) and c) for lab 3

Purpose: to show the conditional numbers for Matricies of size n.
         The solutions are graphed for both the original Differential
         System Matrix as well as the N-1 submatrix

"""

import matplotlib.pyplot as plt
import numpy as np
plt.style.use("ggplot")

def nmatrix(n):
    """
    creates matrix of size n for problem 2
    """
    
    a = np.zeros((n,n+1))
    i,j = np.indices(a.shape)
    a[i==j] = 1.
    a[i==j-1] = -2.
    a[i==j-2] = 1.
    a = np.delete(a,0,-1)
    a = np.array(a)
    return(a)

def smatrix(n):
    """
    creates and array with the boundary conditions 
    substituted into the equations of size n
    """
    
    n = n-1
    a = np.zeros((n,n+2))
    i,j = np.indices(a.shape)
    a[i==j] = 1.
    a[i==j-1] = -2.
    a[i==j-2] = 1.
    top = np.zeros(n+2)
    top[0] = 1.
    bot = np.zeros(n+2)
    bot[n+1] = 1.
    a = np.vstack([top,a,bot])
    
    return(a)

def cond(n,f):
    """
    outputs the conditional number
    for a given function, either,
    nmatrix or fmatrix ("sm", "nm")
    """
    
    if f == "sm":
        A = smatrix(n)
    if f == "nm":
        A = nmatrix(n)
    A = np.linalg.cond(A)
    return(A)
    

# The arrays of values, 
# sequence from 5 to 50
Nval = np.arange(5,51,1)
#  array of conditional numbers for differentail matrix
Ncond  = np.array([cond(i,"nm") for i in Nval])
# array of conditional numbers for sub matrix
SNcond = np.array([cond(i,"sm") for i in Nval])

# Plot it
fig,ax=plt.subplots(1,1,figsize=(5,5))
ax.plot(Nval,Ncond,label="Part (a)")
ax.plot(Nval,SNcond,label="part (b)")
ax.grid(True)
ax.set(xlabel="Maxtrix size 'nxn'",
       ylabel="Conditional Numbers",
       title="Problem 2")
out=ax.legend(loc='upper left')


