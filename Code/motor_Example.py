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
"""

#==============================================================================
# Importing packages
#==============================================================================
import numpy as np
from scipy import sparse
from scipy import linalg
import matplotlib.pyplot as plt
import time as time
from motor import *


#==============================================================================
# Defining motor and simulation parameters
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
mrange = range(2,3)
Frange = np.arange(0.003,0.012,0.003)
Tend = 10

#==============================================================================
# Getting runlength and velocity plots
#==============================================================================

fig1 = getVelVsFload(motor, mrange, Frange, save = False, tol = 1e-6)
fig2 = getRunVsFload(motor, mrange, Frange, Tend, save = False, tol = 1e-6)
