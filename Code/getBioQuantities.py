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
and other functions defined in motor.py

"""

#==============================================================================
# Importing packages
#==============================================================================
import numpy as np
import scipy as sp
from scipy import sparse
from scipy import linalg
from scipy.sparse.linalg import *
import sys
sys.path.append('../')
import copy
from Code.motor import *
import matplotlib.pyplot as plt
import os

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
m = 3
tol = 1e-6
ATP = 0.002
Pback = 2
K_Pd = 0.04 #Original rate from paper
motor = createMotor(l0, Kel, Fs, ATP, Pattach, K_Pd, Pback, motor_type)

# pdf0 = initializePDF(motor,m,Fload,1)
# A,D = getTransitionMatrix(motor, Fload, m, tol, twoSided = False)

#==============================================================================
# Function to obtain average runlengths, cargo loss time, and velocities
#==============================================================================

def get_cargoLossTime(A,pdf0):
    A = A[1::,1::]
    A_T = A.T   
    diag = A_T.diagonal()
    diag_inv = (-1/diag)
    DiagMat = sparse.diags(diag_inv)
    A_T_hat = A_T @ DiagMat
    b = -1*np.ones(A_T_hat.shape[0])  
    A_T_hat = A_T_hat.tocsc()
    cargoLossTime_x, info2 = bicgstab(A_T_hat,b)
    cargoLossTime_x = cargoLossTime_x*diag_inv
    cargoLossTime_x = np.append(0,cargoLossTime_x)
    avg_cargoLossTime = pdf0.dot(cargoLossTime_x)
    
    return cargoLossTime_x, avg_cargoLossTime

def get_cargoLossTime2(A,pdf0):
    # To stabilize this code numerically, we first calculate time to go to 
    # state 1 and then add the time to go from state 1 to state 0. 
    A_T = A.T  
    A_T = A_T[2::,2::]
    diag = A_T.diagonal()
    diag_inv = (-1/diag)
    DiagMat = sparse.diags(diag_inv)
    A_T_hat = A_T @ DiagMat
    b = -1*np.ones(A_T_hat.shape[0])  
    A_T_hat = A_T_hat.tocsc()
    cargoLossTime_x, info2 = bicgstab(A_T_hat,b)
    cargoLossTime_x = cargoLossTime_x*diag_inv
    cargoLossTime_1 = ((cargoLossTime_x @ (A[2::,1]))[0] + 1)/A[0,1]
    cargoLossTime_x = np.append([0,cargoLossTime_1],cargoLossTime_x)
    avg_cargoLossTime = pdf0.dot(cargoLossTime_x)
    
    return cargoLossTime_x, avg_cargoLossTime

def get_runlength(A,D,pdf0):
    A_T = A.T
    D_T = D.T
    b = -1*D_T.sum(axis=1)
    A_T = A_T.tocsc()
    
    runlength_x, info1 = bicgstab(A_T,b)
    avg_runlength = pdf0.dot(runlength_x)
    
    return runlength_x, avg_runlength

def get_runlength2(A,D,pdf0):
    
    A_T = A.T
    D_T = D.T
    A_T = A_T[1::,1::]
    D_T = D_T[1::,1::]
    diag = A_T.diagonal()
    diag_inv = (-1/diag)
    DiagMat = sparse.diags(diag_inv)
    A_T_hat = A_T @ DiagMat
    b = -1*D_T.sum(axis=1)
    A_T_hat = A_T_hat.tocsc()
    runlength_x, info1 = bicgstab(A_T_hat,b)
    runlength_x = runlength_x*diag_inv
    runlength_x = np.append(0,runlength_x)
    avg_runlength = pdf0.dot(runlength_x)
    
    return runlength_x, avg_runlength


def get_runlength3(A,D,pdf0):
    
    # To stabilize this code numerically, we first calculate expected distance  
    # to go to state 1 from any state j and then add the distance to go from 
    # state 1 to state 0. These modified calculations are saved in the pdf in 
    # in bi-directional folder
    A_T = A.T
    D_T = D.T
    
    A_T = A_T[2::,2::]
    D_T = D_T[2::,2::]
    
    diag = A_T.diagonal()
    diag_inv = (-1/diag)
    DiagMat = sparse.diags(diag_inv)
    A_T_hat = A_T @ DiagMat
    b = -1*D_T.sum(axis=1)
    A_T_hat = A_T_hat.tocsc()
    runlength_x, info1 = bicgstab(A_T_hat,b)
    runlength_x = runlength_x*diag_inv
    runlength_1 = ((runlength_x @ (A[2::,1]))[0] + D[1,1])/A[0,1]
    runlength_x = np.append([0,runlength_1],runlength_x)
    avg_runlength = pdf0.dot(runlength_x)
    
    return runlength_x, avg_runlength


def get_SD(A): #Get stationary distribution given rate matrix
#this pdf has 0 appended to it in the beginning
    A3 = A[1::,1::]
    A3.setdiag(np.zeros(A3.shape[0]))
    diag = np.ravel(-1*A3.sum(0)) #Column sum of sparse matrix A
    A3.setdiag(diag) #Setting the diagonal to be sum of column without a for loop
    A3[-1,:] = 1
    b = np.zeros(A3.shape[0])
    b[-1] = 1
    A3 = A3.tocsc()
    pdf, info = bicgstab(A3,b)
    pdf = np.append(0,pdf) # append a zero for the zero motor attached state i.e. code = 1. Just for keeping the size of arrays consitent
    pdf[pdf<0] = 0 # for numerical stability
    
    return pdf

def get_stableSD(A):
    
    A3 = A[1::,1::]
    A3.setdiag(np.zeros(A3.shape[0]))
    diag = np.ravel(-1*A3.sum(0)) #Column sum of sparse matrix A
    A3.setdiag(diag) #Setting the diagonal to be sum of column without a for loop
    # A3[-1,:] = 1
    diag = A3.diagonal()
    diag_inv = (-1/diag)
    D = sparse.diags(diag_inv)
    A4 = A3 @ D
    # A4 = A3 * diag_inv
    A4[-1,:] = 1
    b = np.zeros(A4.shape[0])
    b[-1] = 1
    A4 = A4.tocsc()
    pdf, info = bicgstab(A4,b)
    stablePDF = pdf*(-1/diag)
    stablePDF = stablePDF/sum(stablePDF)
    stablePDF = np.append(0,stablePDF)
    stablePDF[stablePDF<0] = 0
    
    return stablePDF

def get_pi_hat(A):
    
    A_dense = A.todense()
    diag = A.diagonal()
    diag_inv = (-1/diag)
    DiagMat = sparse.diags(diag_inv)
    DiagDense = DiagMat.todense()
    
    R = A @ DiagMat  #calculating routing matrix
    R.setdiag(np.zeros(R.shape[0]))
    
    R[1,0] = 1
    # R[0,0] = 1
    
    I = sparse.identity(R.shape[0])
    R = R - I
    R[-1,:] = 1
    R_dense = R.todense()
    
    b = np.zeros(R.shape[0])
    b[-1] = 1
    
    pi_hat,info_hat = bicgstab(R,b,tol=1e-08)
    pi_hat[pi_hat<0] = 0
    
    N_hat = (pi_hat/pi_hat[0])[1::]
    
    T_hat = N_hat*diag_inv[1::]
    
    T_hat_frac = T_hat/sum(T_hat)
    
    return pi_hat, N_hat, T_hat, T_hat_frac

def get_avgVelocity(A,D):
    
    pdf1 = get_stableSD(A) #get stationary distribution   
    D = D.tocsc()
    avg_vel = np.sum(D,axis = 0).dot(pdf1)
    
    return avg_vel, pdf1


def get_avgMotors(motor,m,Fload,pdf): #under conditional SD that atleast one motor is present because the first entry of PDF is set to 0 atm
    
    n = Nlocations(m, motor, Fload) # Finding Bound on extent, n
    N = Nconfiguration_detach(n,m) #number of rel configurations
    Nm_arr = []
    
    for code in range(1,N+1):
        state, Nm = code2state(code,n) #Here Nm is the number of motors in each state
        Nm_arr.append(Nm)
    
    Nm_arr = np.array(Nm_arr)   
    avg_motors = Nm_arr.dot(pdf)
    
    return avg_motors

#%% Getting Velocities, avg motors, runlengths and cargoLossTime
m=4
K_Pd = 0.04
ATP = 0.002
Pattach = 5
Pback = 2
mrange = [2,3,4]
Fload = 0.00002
# for m in range(3,4):
ATP_range = [0.0002,0.0003,0.0004,0.0005,0.0006,0.0007]
Pback_range = [0.5]
Pback_range = [1,1.5,2,2.5,3,3.5,4,4.5,5,5,5.5,6]
Pback_range = [0.0001,0.001,0.01,0.05,0.1,1,5,10]
# Pback_range = [0.001]
K_Pd_range = [0.0005,0.001,0.005,0.01,0.05,0.1]
# Pattach_range = [1,10,20,30,40,50]
# Pback_range = 0.1*np.array([1,5,10,15])
# Pback_range = np.array([0.5,1,1.5,2,2.5,3,3.5,4,4.5,5,5.5,6,6.5,7,7.5,8,8.5,9,9.5,10])
# ATP = 0.0005
# ATP_range = [0.5e-4,0.6e-4,0.7e-4,0.8e-4,0.81e-4,0.84e-4,0.9e-4,1e-4,2e-4,10e-4,50e-4,100e-4]
dF = 0.002
ATP_range = [0.5e-4,0.75e-4,1e-4,2e-4,10e-4,50e-4]
Frange = np.arange(0.002, 0.022, dF)
motor1 = createMotor(l0,Kel,Fs,ATP,Pattach,K_Pd,Pback,motor_type)

for m in range(1,4):    
    start1 = time.time()
    avg_vel_1sided = []
    avg_vel_2sided = []
    avg_runlength_1sided = []
    avg_runlength_2sided = []
    avg_cargoLossTime_1sided = []
    avg_cargoLossTime_2sided = []
    PDF_1sided_arr = []
    PDF_2sided_arr = []
    avgNm_1sided_arr = [] 
    avgNm_2sided_arr = [] 
    
    
    print(m)

    for Fload in Frange:
        start2 = time.time()
        print(Fload)
        
        
        A1, D1 = getTransitionMatrix(motor1, Fload, m, tol, twoSided = False)
        # A2, D2 = getTransitionMatrix_relSpace(motor1, Fload, m, bigRelSpace_m4, tol, twoSided = True)
        pdf0 = initializePDF(motor1,m,Fload,1)
        # pdf1 = np.append(0,get_pi_hat(A1)[3]) #using modified routing matrix method to get fraction of time spent in each state
        # pdf2 = np.append(0,get_pi_hat(A2)[3])
        
        
        
        # avg_vel1, pdf1 = get_avgVelocity(A1,D1)
        # avg_vel2, pdf2 = get_avgVelocity(A2,D2)
        
        # PDF_1sided_arr.append(pdf1)
        # PDF_2sided_arr.append(pdf2)
        
        # avg_vel_1sided.append(avg_vel1)
        # avg_vel_2sided.append(avg_vel2)
        
        # avgNm_1sided_arr.append(get_avgMotors(motor,m,Fload,pdf1))
        # avgNm_2sided_arr.append(get_avgMotors(motor,m,Fload,pdf2))
        
        runx_x_1, avg_run1 = get_runlength3(A1,D1,pdf0)
        # runx_x_2, avg_run2 = get_runlength3(A2,D2,pdf0)
        
        avg_runlength_1sided.append(avg_run1)
        # avg_runlength_2sided.append(avg_run2)
        
        hitTime_x_1, avg_time1 = get_cargoLossTime2(A1,pdf0)
        # hitTime_x_2, avg_time2 = get_cargoLossTime2(A2,pdf0)
        
        avg_cargoLossTime_1sided.append(avg_time1)
        # avg_cargoLossTime_2sided.append(avg_time2)
        
        avg_vel1 = avg_run1/avg_time1
        # avg_vel2 = avg_run2/avg_time2
        
        avg_vel_1sided.append(avg_vel1)
        # avg_vel_2sided.append(avg_vel2)
        
        
    avg_vel_1sided = np.array(avg_vel_1sided).ravel()
    # avg_vel_2sided = np.array(avg_vel_2sided).ravel()
    avg_runlength_1sided = np.array(avg_runlength_1sided).ravel()
    # avg_runlength_2sided = np.array(avg_runlength_2sided).ravel()
    avg_cargoLossTime_1sided = np.array(avg_cargoLossTime_1sided).ravel()
    # avg_cargoLossTime_2sided = np.array(avg_cargoLossTime_2sided).ravel()
        
    
    end1 = time.time()  
    print('time elapsed for m=', m, 'is', "%.2f" % (end1-start1), 'seconds')
    
    #=========================================================================
    # Plotting Routines
    plt.rc('xtick', labelsize=20) 
    plt.rc('ytick', labelsize=20) 
    
    
    xvar = Frange
    xlabel = ' Fload (pN) '
    legendVar1 = m
    legendLabel1 = ' motors '
    legendVar2 = K_Pd
    legendLabel2 = ' K_Pd '
    legendVar3 = ATP
    legendLabel3 = ' ATP '
    legendVar4 = Pback
    legendLabel4 = ' Pback '
    legend = legendLabel1 +'= {}'.format(str(legendVar1)) 
    subtitle = legendLabel2 +'= {}'.format(str(legendVar2) + legendLabel3 +'= {}'.format(str(legendVar3) + legendLabel4 +'= {}'.format(str(legendVar4))))
    
    fntSz1 = 20
    fntSz2 = 16
    
    
    # fig1,ax1 = plt.subplots(1)
    # ax1.plot(xvar, avg_vel_1sided, color="blue", label = '1 sided ' + legend)
    # ax1.plot(xvar, avg_vel_2sided, color="orange", label = '2 sided ' + legend)
    # plt.grid(True)
    # plt.title(subtitle, fontsize = fntSz2)
    # fig1.suptitle('Velocity vs ' + xlabel , fontsize = fntSz1)
    # ax1.set_xlabel(xlabel, fontsize = fntSz1)
    # ax1.set_ylabel('Average Velocity (nm/s)', fontsize = fntSz1)
    # ax1.legend(fontsize = fntSz1)    
    
    # ax2 = ax1.twinx()
    # ax2.scatter(xvar,avgNm_1sided_arr,color="green",marker="o",label = '1 sided motors')
    # ax2.scatter(xvar,avgNm_2sided_arr,color="black",marker="o", label = '2 sided motors')
    # ax2.set_ylabel('Average Number of Motors',fontsize = fntSz1)
    # ax2.legend(fontsize = fntSz1) 
    
    fig1 = plt.figure(1)
    plt.plot(xvar, avg_vel_1sided, marker="*", label = '1 sided ' + legend)
    plt.grid(True)
    plt.title(subtitle, fontsize = fntSz2)
    fig1.suptitle('Velocity vs ' + xlabel, fontsize = fntSz1)
    plt.xlabel(xlabel, fontsize = fntSz1)
    plt.ylabel('Velocity (nm/s)', fontsize = fntSz1)
    plt.legend( fontsize = fntSz1)
    
    fig2 = plt.figure(2)
    plt.plot(xvar, avg_runlength_1sided, marker="*", label = '1 sided ' + legend)
    # plt.plot(xvar, avg_runlength_2sided, color="orange", label = '2 sided ' + legend)
    plt.grid(True)
    plt.title(subtitle, fontsize = fntSz2)
    fig2.suptitle('Runlength vs ' + xlabel, fontsize = fntSz1)
    plt.xlabel(xlabel, fontsize = fntSz1)
    plt.ylabel('Average runlength (nm)', fontsize = fntSz1)
    plt.legend( fontsize = fntSz1)   
    
    fig3 = plt.figure(3)
    plt.plot(xvar, avg_cargoLossTime_1sided, marker="*", label = '1 sided '+ legend)
    # plt.plot(xvar, avg_cargoLossTime_2sided, label = '2 sided '+ legend)
    plt.grid(True)
    plt.title(subtitle, fontsize = fntSz2)
    fig3.suptitle('Cargo Loss Time vs ' + xlabel, fontsize = fntSz1)
    plt.xlabel(xlabel, fontsize = fntSz1)
    plt.ylabel('Average cargo Loss time (s)', fontsize = fntSz1)
    plt.legend( fontsize = fntSz1)   
    
