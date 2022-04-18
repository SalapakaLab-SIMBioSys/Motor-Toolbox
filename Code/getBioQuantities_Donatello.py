# -*- coding: utf-8 -*-
"""
Created on Mon Jan 27 14:59:54 2020

Author: Rachit Shrivastava

MEASURE SYSTEM
nm for length
M for concentration
Kg for weight
s for second
K for temperature
nN for force

This file contains functions to extract biologically relevant quantities
from single motor properties. This imports motor.py and operates on motor class
and other functions defined in motor.py. Exactly as coded by Donatello in Matlab. Results match for both vel and runlength.

"""

#==============================================================================
# Importing packages
#==============================================================================
import numpy as np
import scipy as sp
from scipy import sparse
from scipy import linalg
import copy
from motor_Donatello import *
# from motor_modified_copy import *
import matplotlib.pyplot as plt

#import pylab as pl
import time as time

#==============================================================================
# Initializing motor
#==============================================================================
l0 = 110
Kel = 0.32e-3
ds = 8
Fs = 0.006
Pattach = 5
motor_type = 'kinesin'
Fload = 0.002
m = 2
tol = 1e-6
ATP = 0.002
K_Pd = 0.04
Pback = 2
motor = createMotor(l0, Kel, Fs, ATP, Pattach, motor_type)
# motor = createMotor(l0,Kel,Fs,ATP,Pattach,K_Pd,Pback,motor_type)


#%%
#==============================================================================
# Function to obtain average runlengths and its variation as load force is
# varied
#==============================================================================

dF = 0.002
dt = 0.002

#Frange = np.arange(0.0002, m*motor.Fs - 0.0002, dF)
Frange = np.arange(0.002, 0.012, dF)
# Frange = [0.002]
trange = np.arange(0,100,dt)
#trange = np.arange(0.003,0.006,dt)

for m in range(2,4):

    start = time.time()
    runLength_Tinf = np.empty(0)
        
    for Fload in Frange:
        A, D = getTransitionMatrix(motor, Fload, m, tol)
        # A,D = getTransitionMatrix_1sided(motor, Fload, m, tol)
        pdf0 = initializePDF(motor,m,Fload,1)
        pdf = pdf0
        x = 0
        A = A.tocsr()
        D = D.tocsr()
        
        for t in trange:
            
        #        x += np.sum(D,axis=0) @ pdf * dt
        #        x += np.dot(np.sum(D,axis=0) , pdf) * dt
            x += np.sum(D,axis=0).dot(pdf) * dt
        #        pdf += A @ pdf * dt
            pdf += A.dot(pdf) * dt
            
            
        #        
               
        runLength_Tinf = np.append(runLength_Tinf, x)
    
    end = time.time()  
    print('time elapsed for motor number =', m, 'is', "%.2f" % (end-start), 'seconds')
    
    fig = plt.figure(91)
    plt.plot(Frange, runLength_Tinf, marker = "*", label = 'No. of Motors = {}'.format(str(m)))
    plt.grid(True)
    fig.suptitle('Runlength vs Load Force')
    plt.xlabel('Force (pN)')
    plt.ylabel('Runlength (nm)')
    plt.legend()    
    
 
#%% Getting Velocities

for m in range(2,4):
    start = time.time()
    avg_vel = np.empty(0)
    print(m)

    dF = 0.002
    ATP_range = [0.5e-4,0.75e-4,1e-4,2e-4,10e-4,50e-4]
    Frange = np.arange(0.002, 0.012, dF)
    
    
    for Fload in Frange:
       
        A, D = getTransitionMatrix(motor, Fload, m, tol)
        # A,D = getTransitionMatrix_1sided(motor, Fload, m, tol)
        A2 = A[1::,1::]
        
        for ii in range(0,np.size(A2[0])):
            A2[ii,ii] = 0
            
        for ii in range(0,np.size(A2[0])):
            A2[ii,ii] = -np.sum(A2[:,ii])
            
        A2 = A2.todense()
        
        pdf = linalg.null_space(A2)
        pdf = pdf/np.sum(pdf)
        pdf = np.append(0,pdf)
        
        D = D.todense()
        
        avg_vel = np.append(avg_vel, np.sum(D,axis = 0).dot(pdf))
        
    end = time.time()  
    print('time elapsed for motor number =', m, 'is', "%.2f" % (end-start), 'seconds')
    
    fig = plt.figure(92)
    plt.plot(Frange, avg_vel, marker = "*", label = 'No. of Motors = {}'.format(str(m)))
    plt.grid(True)
    fig.suptitle('Velocity vs Load Force')
    plt.xlabel('Force (pN)')
    plt.ylabel('Average Velocity (nm)')
    plt.legend()    
    