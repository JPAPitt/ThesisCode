"""
Created on Wed Apr 22 10:09:21 2015

@author: jordan

Python file that allows the use of the C library (LaxWendroffSolver.c) to solve the Serre equuations
"""

from LaxWendroffSolver import *
from scipy import *
import csv
import os
from pylab import plot, show, legend,xlim,ylim,title,xlabel,ylabel


#Python-C interface functions

def copyarraytoC(a):
	#Return a copy of Python list (a) as a C list
    n = len(a)
    b = mallocPy(n)
    for i in range(n):
        writetomem(b,i,a[i])
    return b
  
def copyarrayfromC(a,n):
	#Return a copy of C list 'a' with length n as a Python list 
    b = [0]*n
    for i in range(n):
        b[i] = readfrommem(a,i)
        
    return b

def makevar(sx,ex,dx,st,et,dt): 
	#Make the x and t lists
	#Given the start location sx and end location ex and resolution dx 
	#Given the start time st and end time et and resolution dt 
    x = arange(sx, ex, dx)
    t = arange(st, et, dt)
    
    return x,t 
        

def sech2 (x):
  #returns sech^2(x)
  a = 2./(exp(x) + exp(-x))
  return a*a

def soliton (x,t,g,a0,a1):
  #returns h(x,t) for the soliton example with water depth a_0, amplitude a_1 and gravitational acceleration g
  c = sqrt(g*(a0 + a1))
  phi = x - c*t;
  k = sqrt(3.0*a1) / (2.0*a0 *sqrt(a0 + a1))
  return a0 + a1*sech2(k*phi)
  
def solitoninit(n,a0,a1,g,x,t0,dx):
    #returns vectors h(x,t), u(x,t) for the soliton example with water depth a_0, amplitude a_1 and gravitational acceleration g
    
    h = zeros(n)
    u = zeros(n)
    c = sqrt(g*(a0 + a1))
    for i in range(n):
        h[i] = soliton(x[i],t0,g,a0,a1)
        u[i] =  c* ((h[i] - a0) / h[i])
    
    return h,u


#---------- Soliton Example ----------------

#physical quantities
g = 9.81
a0 = 1.0
a1 = 0.7


#temporal and spatial resolution
dx = 100.0 / (2**12)

#CFL condition
Cr = 0.5
l = 1.0 / (sqrt(g*(a0 + a1)))
dt = Cr*l*dx


#start and end times and locations
startx = -50.0
endx = 250.0 + dx
startt = 0
endt = 50 + dt


#generate spatial and temporal space    
x,t = makevar(startx,endx,dx,startt,endt,dt)
n = len(x)

# vectors h and u at start time, and ph and pu at start time - dt
h,u = solitoninit(n,a0,a1,g,x,0.0,dx)
ph,pu = solitoninit(n,a0,a1,g,x,-dt,dx)


# Boundary conditions (note that nBCs > nBC)   
nBC = 3
nBCs = 4
u0 = zeros(nBCs)
u1 = zeros(nBCs)    
h0 = a0*ones(nBCs)
h1 = a0*ones(nBCs)
    
#make all the required vectors, C vectors
h_c = copyarraytoC(h)
u_c = copyarraytoC(u)
pubc_c = copyarraytoC(concatenate([u0[-nBC:],pu,u1[:nBC]]))
phbc_c = copyarraytoC(concatenate([h0[-nBC:],ph,h1[:nBC]]))
h0_c  = copyarraytoC(h0)
h1_c  = copyarraytoC(h1)
u0_c  = copyarraytoC(u0)
u1_c  = copyarraytoC(u1)

      
for i in range(1,len(t)):            
    evolvewrap(u_c, h_c, pubc_c, h0_c, h1_c,u0_c, u1_c,g,dx,dt,nBC, n,nBCs)    
    print (t[i])

#Get h and u from C lists to python list
uFinal = copyarrayfromC(u_c,n)
hFinal = copyarrayfromC(h_c,n)


#get the exact solutions for h and u at the final time
hAnalytic,uAnalytic = solitoninit(n,a0,a1,g,x,t[i],dx)

#Plot h, both analytic and numerical
plot(x,hAnalytic,label='Analytic')
plot(x,hFinal,label='Numerical')
title('Soliton Problem Example')
xlabel('x (m)')
ylabel('h (m)')
legend()

#Free all the allocated memory
deallocPy(u_c)   
deallocPy(h_c)
deallocPy(h0_c)
deallocPy(h1_c)
deallocPy(u0_c)
deallocPy(u1_c)



