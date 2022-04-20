#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Code for running the Motor-Toolbox and generating runlength vs Fload and 
velocity vs Fload plots

MEASURE SYSTEM
nm for length
M for concentration
Kg for weight
s for second
K for temperature
nN for force

Define the motor parameters:
    
    l0  - type: float - rest length of motor in nm
    
    Kel - type: float - elasticity coefficient of motor in nN/nm
    
    ds  - type: int - step size of motor in nm
    
    Fs - type: float - stall force of motor in nN
    
    Pattach -type: int - reattachment rate of motor in per second
    
    Pback - type: float - detachment rate of motor after stall in per second
    
    K_Pd - type: float - constant with which stepping rate is multiplied to get detachment rate
    
    motor_type - type: string - motor type, 'kinsin', 'myosin', 'dynein'
    
Define environment parameters:
    
    ATP - type: float - ATP concentration in M
    

Define Simulation Parameters:    

    motor - object of class motorClass - example: motor = createMotor(l0, Kel, Fs, ATP, Pattach, motor_type)
    
    mrange - type: array of integeres - example: mrange = range(2,4)
    
    Frange - type: numpy array - example: Frange = np.arange(0.002,0.012,0.002) 
        - Range of load force user wants to simulate, must start with a non zero float
    
    save - type: Boolian - User input to save figure generated
    
    Tend - type: int - Simulation time horizon
   
    tol - tolerance for calculating equilibrium position of cargo, default = 1e-6
    
   

"""

#==============================================================================
# Importing packages
#==============================================================================
import numpy as np
from motor import *

#==============================================================================
# Defining motor, environment, and simulation parameters
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
mrange = range(2,4)
Frange = np.arange(0.002,0.012,0.0005)
Tend = 10

#==============================================================================
# Getting runlength and velocity plots
#==============================================================================

fig1 = getVelVsFload(motor, mrange, Frange, save = False, tol = 1e-6)
fig2 = getRunVsFload(motor, mrange, Frange, Tend, save = False, tol = 1e-6)
