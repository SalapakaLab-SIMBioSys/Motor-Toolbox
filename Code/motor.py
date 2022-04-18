# -*- coding: utf-8 -*-
"""
Created on Mon Dec 30 16:23:21 2019

Author: Rachit Shrivastava

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
import scipy as sp
from scipy import sparse
import copy

#import pylab as pl
#import matplotlib.pyplot as plt
#import time as time


#==============================================================================
# Class for motors
#==============================================================================

class motorClass:
    #==========================================================================
    # Class Variables
    #==========================================================================
    
    # Physical parameters
    Kb = 1.3806503e-5; # Boltzmann constant in nm^2 kg s^(-2) K^(-1)
    T = 300 # Temperature in Kelvin
    # ATP = 2e-3 # ATP concentration in M

    # Bio-chemo-mechanical motion parameters
    kcat = 105 # s^(-1)
    kon = 2e6 # M^(-1) s^(-1)
    k0off = 55 # s^(-1)
    
    # Microtubule properties and binding Bio-chemo-mechanical parameters
    ds = 8 # Periodicity of MT filament in nm
    dl = 1.6 # nm find explanation for this
    delta_l = 1.3 # nm find explanation for this
    A = 107 # Mean step number
    B = 0.029e-6 # M
    Pback = 2 # Find explanation for this
    
    # Bead parameters (May be used later)
    bead_r = 250 # Bead radius in nm
    eta = 0.89e-12; # Coefficient of viscosity for water in Kg/(s nm)
    Xi = 6*np.pi*eta*bead_r # Friction force Kg/s
    D = Kb*T/Xi ; # nm^2/s
    
    #==========================================================================
    # Init Function
    #==========================================================================
    
    def __init__(self, l0, Kel, Fs, ATP, Pattach, motor_type):
        
        self.l0 = l0
        self.Kel = Kel
        self.Fs = Fs
        self.ATP = ATP
        self.Pattach = Pattach
        self.motor_type = motor_type
   
    #==========================================================================
    # Probability of stepping
    # (this is dependent on the force the motor is experiencing)
    #==========================================================================
    
    def Pstep(self, Fm):
        if Fm > self.Fs: # Force exerted due to stretching the linkage > Stalling force Fs 
            P = 0 # Probability of stepping is 0 for force beyond stall force
            epsilon = 0
            return [P, epsilon]
        
        if Fm < 0: # Force exerted due to linkage stretching < 0 --> it is taken equal to 0 i.e. "not stretched" stage
            Fm = 0
        
        koff = self.k0off * np.exp( (Fm*self.dl)/(self.Kb*self.T) ) # s^(-1) ; As per formula of Koff, from ref.paper[7]
        km = (self.kcat+koff)/self.kon # M
        epsilon = 1-(Fm/self.Fs)**2 # dimensionless ; epsilon = epsilon(F) i.e. prob. of binding with MT   
        
        P = (self.kcat*self.ATP)/(self.ATP+km); # s^(-1)
        
        if Fm > 0: # i.e. Force exerted due to linkage stretching is 0<F<Fs
            P = P*epsilon
        
        return [P, epsilon]
    
    #==========================================================================
    # Probability of detachment 
    # (this is dependent on the force the motor is experiencing)
    #==========================================================================
    
    def Pdetach(self, Fm):
        if abs(Fm) >= self.Fs:
            P = self.Pback
        else:
            #P1 = self.B/self.ATP #This is not being used anythere.
            P2 = (1/self.A)* np.exp(  (0.008/self.Fs)*(abs(Fm)*self.delta_l)/(self.Kb*self.T) )
            
            P = (((self.ATP/P2)/(self.ATP +self.B*(1+1/P2)))**(-1)) * self.Pstep(Fm)[0]
        return P
      
    def Pdetach_new(self, Fm):
        if abs(Fm) >= self.Fs:
            P = self.Pback
        else:
            P = 0.04 * self.Pstep(Fm)[0]
            
        return P
    
    #==========================================================================
    # Efficiency (this is dependent on the force the motor is experiencing)
    #==========================================================================
    
    def efficiency(self, Fm):
        if Fm < 0:
            Fm = 0
            
        e = 1-(Fm/self.Fs)**2
        
        
        
        if Fm > self.Fs:
            e = 0
        
        e = e*90/100+10/100*np.exp(-2*(Fm/self.Fs)**2) # Comment if not required. Use if a more general expression for efficiency is needed.
        
        return e
    
    #==========================================================================
    # Processivity (this is dependent on the force the motor is experiencing)
    #==========================================================================
    def processivity(self, Fm):

        P2 = (1/self.A)*np.exp( (0.008/self.Fs)*(abs(Fm)*self.delta_l)/(self.Kb*self.T) )
        L = (self.ds*self.ATP/P2)/(self.ATP +self.B*(1+1/P2))
        
        return L
    
    #==========================================================================
    # Spring Force (this is dependent on how much the motor is stretched)
    #==========================================================================
    def springForce(self, l):

        F = 0
        if l > self.l0:
            F = self.Kel * (l - self.l0)
        
        if l < -1*self.l0:
            F = self.Kel * (l + self.l0)
            
        return F
    
#=============== End of Class motorClass ======================================
        
    
#==============================================================================
# Defining other relevant functions here
#==============================================================================

#==============================================================================
# Function to get a column from a list of lists. This comes in handy when 
# working with lists instead of matrix
#==============================================================================
def getColumn(lol,colNum):
    column = []
    for ii in range(0,5):
        column.append(lol[ii][colNum])
        
    return column

#==============================================================================
# Function to create an object of motorClass
#==============================================================================
def createMotor(l0, Kel, Fs, ATP, Pattach, motor_type):
   k = motorClass(l0, Kel, Fs, ATP, Pattach, motor_type) 
   return k

#==============================================================================
# Function to find number of locations (or extent) given max number of motors,
# rest length, periodicity of filament, stall force, and spring constant
#==============================================================================
def Nlocations_old(m,l0,ds,Fs,Kel,Fm):
    
    dmax1 = (m*Fs-Fm)/Kel + ds
    dmax2 = Fm/Kel
    dmax = max(dmax1,dmax2) + 2*l0

    n = np.ceil(dmax/ds) + 1
    n = n + 1 # Figure out +1 thing here.
    
    return int(n)

def Nlocations(m, motor, Fm):
    
    dmax1 = (m*motor.Fs-Fm)/motor.Kel + motor.ds
    dmax2 = Fm/motor.Kel
    dmax = max(dmax1,dmax2) + 2*motor.l0

    n = np.ceil(dmax/motor.ds) + 1
    n = n + 1 # Figure out +1 thing here.
    
    return int(n)

#==============================================================================
# Function to determine number of configurations EXACTLY given m motors 
# and n locations
#==============================================================================
def Nconfiguration(n,m):
    if m < 0:
        N = 0
        
    elif m == 0:
        N = 1
        
    else:
        N = 1
        for ii in range (n,n+m-1):
            N = N*ii
        N = N/ np.math.factorial(m-1)
    return int(N)

    #Easier Implementation
    # N = np.math.factorial(n+m-2)/(np.math.factorial(n-1)*np.math.factorial(m-1))
    
#==============================================================================    
# Function to determine total number of relative configurations EXACTLY
# given atmost m motors on cargo and n locations
#==============================================================================
def Nconfiguration_detach(n,m):
    N = 0
    for ii in range (0,m+1):
        N = N + Nconfiguration(n,ii)
        
    return int(N)


#==============================================================================    
# Function to determine the exact state given the code of the state
# Inputs are code = value from 1 to N (for total N relative configurations)
# i.e. code = no. of some relative configuration
# and n = bound on extent
#==============================================================================
def code2state(code,n):
    
    Ns = n # Bound on extent
    m = 0
    
    while Nconfiguration_detach(n,m) < code:
        m += 1
        
    Nm = m
    state = np.zeros(Ns, dtype = 'int')
    
    if Nm == 0:
        return state, Nm
    
    
    if Nm != 0:
        
        code = code - Nconfiguration_detach(n,m-1) - 1
        n = n-1
        
        k = 1
        state[0] = 1
        m = m-1
        
        if m > 0:
            while k <= Ns:
                N = Nconfiguration(n+1,m)
                if code >= N:
                    k = k+1
                    n = n-1
                    code = code -N
                else:
                    m = m-1
                    state[k-1] = state[k-1] + 1
                    
                if m <= 0:
                    break
        state = np.array(state, dtype='int')
        return state, Nm


#==============================================================================
# Trim State function to obtain a trimmed state
#==============================================================================   
def trimState(state):
    if np.sum(state) == 0:
        return state
    
    else:
        while state[0] == 0:
            state = np.append(state[1:], 0)
        
        state = np.array(state, dtype='int')    
        return state


#==============================================================================    
# Function to determine the code for a state given the state
# Input is a state vector
        
# N = enumerate_combination_rep(s)
# s is a vector of non-negative integers representing the number of
# agents in Ns positions; Ns=n from before
# m represents the maximum number of agents
# a combination of (#'0's) objects chosen with repetition from (#'1's+1) positions
#==============================================================================
def state2code(state):
    Ns = len(state)
    n = Ns
    m = np.sum(state)
    code = Nconfiguration_detach(n,m-1)
    
    state = trimState(state)
    
    m = m-state[0]
    n = n-1
    
    for k in range(2,Ns+1):
        code = code + (m>0)*Nconfiguration(n+1,m)
        m = m - state[k-1]
        n = n-1
        
    code = code + 1
    
    return int(code)

#==============================================================================    
# Function that gives the positions of motors on microtubule given a state
#==============================================================================
def state2pos(state):
    pos = np.empty(0)

    for ii in range (1,len(state)+1):
        if state[ii-1] > 0 :
            pos = np.append(pos, ii * np.ones(int(state[ii-1])))
    
    pos = np.array(pos, dtype='int')    
    return pos

  
#==============================================================================    
# A wrong springForce implementation is used in previous code, this function is
# just recreation of the "wrong" function in python in order to check the 
# validity of other functions like eqm_cargo and state2pos
#==============================================================================
def dummySpringForce(l,Kel,l0):

        F = 0
        if l > l0:
            F = Kel * (l - l0)
        
        if l < -l0:
            F = Kel * (l + l0)
            
        return F

#==============================================================================
# Approx Function for determining cargo equilibrium position given the  
# configuration of motors with no load force
# This is the corrected version of eqm_cargo, however it is not used. 
#==============================================================================     
def eqm_cargo(motor, p, Fload, tol = 1e-6):
    p = p*motor.ds
    p = np.sort(p)
    n = len(p)
    
    x0 = (Fload + motor.Kel*n*motor.l0 + motor.Kel*(np.sum(p)))/n #First Guess
    
    #Check net force:
    nf = Fload
    x = x0
    for ii in range(0,n):
        nf = nf + motor.springForce(x0 - p[ii])
        
    if nf >= 0:
        xhigh = x
        xlow = xhigh - nf/motor.Kel
        
    else:
        xlow = x
        xhigh = xlow + abs(nf)/motor.Kel
        
    while abs(nf) > tol:
        
        x = (xhigh + xlow)/2
        
        nf = Fload
        
        for jj in range(0,n):
            nf = nf + motor.springForce(x -p[jj])
            
        if nf > 0:
            xhigh = x
            
        else:
            xlow = x
            
    return x

#==============================================================================
# Approximate Function for determining cargo equilibrium position given  
# the configuration of motors with no load force
# This is the"wrong" version of eqm_cargo which is the same as 
# Donatello's implementation. This uses dummySpringForce function and this was 
# written just to check the validity of python engine.
#============================================================================== 
def dummy_eqm_cargo(motor, p, Fload, tol = 1e-6):
    p = motor.ds*p
    p = np.sort(p)
    n = len(p)
    
    x0 = (Fload + motor.Kel*n*motor.l0 + motor.Kel*(np.sum(p)))/n #First Guess

    #Check net force:
    nf = Fload
    x = x0
    for ii in range(0,n):
        nf = nf + dummySpringForce(x0 - p[ii],motor.l0,motor.Kel)
        
    if nf >= 0:
        xhigh = x
        xlow = xhigh - nf/motor.Kel
        
    else:
        xlow = x
        xhigh = xlow + abs(nf)/motor.Kel
        
    while abs(nf) > tol:
        
        x = (xhigh + xlow)/2
        
        nf = Fload
        
        for jj in range(0,n):
            nf = nf + dummySpringForce((x -p[jj]),motor.l0,motor.Kel)
            
        if nf > 0:
            xhigh = x
         
        else:
            xlow = x
    return x

#==============================================================================
# Exact Function for determining cargo equilibrium position given the  
# configuration of motors. This is the actual function which is used in
# calculating transition matrix i.e. in upcoming function trans_mat_general
    
'''
 synopsis:
 eq=equilibrium_cargo(p,Kel,l0,F0)
 ---
 p is a vector containing the positions of the motors
 Kel is the elastic constant of the linkages
 l0 is the rest length of the linkages
 FLoad is the load applied on the cargo

'''
#============================================================================== 
def equilibrium_cargo(motor, p, Fload, tol = 1e-6):
    p = p*motor.ds
    p = np.sort(p)
    
    N = len(p)
    
    if Fload == 0: 
        #Zero load force on cargo, 'eqm_cargo' or 'dummy_eqm_cargo' only used here
        eq = eqm_cargo(motor, p, Fload, tol)
        return eq
    
    else:
        if Fload > 0: # Positive load force on cargo
            eq = (p[-1]- motor.l0) - Fload/motor.Kel
            imin = 1
            imax = N-1
            
        if Fload < 0: # Negative load force on cargo
            eq = (p[0]+ motor.l0) - Fload/motor.Kel
            imin = 2
            imax = N
        
        while imin <= imax:
            F = 0
            
            for ii in range(0,N):
                if (p[ii] - eq) > motor.l0:
                    F = F + motor.Kel * (p[ii] - eq - motor.l0)                   
                    
                if (p[ii] - eq) < -motor.l0:
                    F = F + motor.Kel * (p[ii] - eq + motor.l0)
                                
            if F < Fload:
                imin = imin + 1
                
            if F > Fload:
                imax = imax - 1
                
            if abs(F-Fload) <= abs(Fload)*tol:
                break
            
            eq = (1/(len(p)-(imax-imin+1)))*(sum(p[:imin-1] + motor.l0) + sum(p[imax:] - motor.l0) - Fload/motor.Kel)

            # check if there are no strecthed connections
#            print(sum(abs(p[imin-1:imax] - eq) > motor.l0))
            if (sum(abs(p[imin-1:imax] - eq) > motor.l0)  == 0):
                break
        
        return eq

#==============================================================================
# Function to calculate transition matrix. This is the main function 
# [A n D Psingle self_state StepStat, ProbCargoStep] = trans_mat_general(motor, F, m)
'''
This function takes in the motor parameters, load force and number of motors
to give out a N by N transition matrix where N is the maximum number of relative
configurations given maximum m number of motors are present on the cargo.

'''
#============================================================================== 
def getTransitionMatrix(motor, Fload, m, tol):
    Psingle = []
    selfState = []
    
    # Finding Bound on extent, n
    n = Nlocations(m, motor, Fload)
    
    # Finding number of possible relative configuration
    N = Nconfiguration_detach(n,m)
    
    A = sparse.dok_matrix((N, N), dtype=np.float32) #Create NxN sparse matrix
    #A = A.tocoo() #Converting from Dictionary of keys to coordinate format matrix
    D = sparse.dok_matrix((N, N), dtype=np.float32) #Create anther NxN sparse matrix
    #D = D.tocoo() 
    stepStat = []
    probCargoStep = []
    
    for code in range(1,N+1): 
        # step through each code ascribed to all the relative configurations
        state, Nm = code2state(code,n) 
        # state corresponding to the current code obtained
        pos = state2pos(state) 
        # positions of motors corresponding to state obtained
        
        if pos.size == 0:
            x_c = 0        
        else:
            x_c = equilibrium_cargo(motor, pos, Fload, tol) 
            #Calculate the equilibrium position of cargo       
        
    #==================== Compute Step Transition rates =======================
        for ii in range(0,n-1):
            if state[ii] > 0:
                Fm = motor.springForce((ii+1)*motor.ds - x_c) 
                #Calculate the force acting on ONE motor at a position
                
                '''
                Sends l= (i*ds-x_c) which is stretched linkage adjusted by 
                cargo position x_c. 
                Fm(i) contains value of force exerted by a motor at i th
                location
                
                '''
                #Step transition happening
                newState = np.copy(state)
                newState[ii] = newState[ii] -1  
                newState[ii+1] = newState[ii+1] + 1 
                #new state obtained by this point
    
                newCode = state2code(newState) 
                #code for new state obtained
                newPos = state2pos(newState) 
                #new motor positions in the new state obtained
                new_x_c = equilibrium_cargo(motor, newPos, Fload, tol) 
                #new equilibrium position obtained
                
                if newCode != code:
                    A[newCode-1,code-1] += state[ii] * (motor.Pstep(Fm)[0]) 
                    #[0] index in Pstep gives the probability and [1] index gives the epsilon
                else:
                    Psingle.append(motor.Pstep(Fm)[0])
                    selfState.append(code)
                    
                #Not used immediately in pdf calculation (?)
                D[newCode-1,code-1] += (new_x_c - x_c) * state[ii] * (motor.Pstep(Fm)[0])
                tempVar = np.array([code,newCode,1])
                if new_x_c != x_c:
                    stepStat.append(tempVar)
                    probCargoStep.append([(new_x_c - x_c), (motor.Pstep(Fm)[0]) * state[ii]])
              
    #==================== Compute Detachmemt Transition rates =================
        for ii in range(0,n):
            if state[ii] > 0:
                Fm = motor.springForce((ii+1)*ds - x_c)
                newState = np.copy(state)
                newState[ii] = newState[ii] - 1 
                # detachment event occured, new state obtained
                newCode = state2code(newState) #code for new state obtained
                newPos = state2pos(newState)
                
                if newPos.size == 0:
                    new_x_c = x_c
                else:
                    new_x_c = equilibrium_cargo(motor, newPos, Fload, tol)
                
                if newCode != code:
                    A[newCode-1,code-1] += state[ii] * motor.Pdetach(Fm)
                
                #Not used immediately in pdf calculation (?)
                D[newCode-1,code-1] += (new_x_c - x_c) * state[ii] * motor.Pdetach(Fm)
                tempVar = np.array([code,newCode,-1])
                if new_x_c != x_c:
                    stepStat.append(tempVar)
                    probCargoStep.append([(new_x_c - x_c),  state[ii] * motor.Pdetach(Fm)])  
                
    #==================== Compute attachment Transition rates =====================
        floatMotors = int(m - sum(state))
        
        if sum(state) < m and sum(state) > 0:
            minReattachLoc = int(np.ceil((x_c - motor.l0)/motor.ds))
            maxReattachLoc = int(np.floor((x_c + motor.l0)/motor.ds))
            
            if minReattachLoc <= 0:
                shiftedState = np.append(np.zeros(1-minReattachLoc) ,state)
                shiftedState = shiftedState[:n]
                maxReattachLoc += 1 - minReattachLoc
                minReattachLoc = 1    
            else:
                shiftedState = np.copy(state)
                
            for jj in range(minReattachLoc-1, maxReattachLoc):
                newState = np.copy(shiftedState)
                newState = np.array(newState, dtype = 'int')
                newState[jj] += 1
                
                shiftedPos = state2pos(shiftedState)
                shifted_x_c = equilibrium_cargo(motor, shiftedPos, Fload, tol)
                
                newPos = state2pos(newState)
                new_x_c = equilibrium_cargo(motor, newPos, Fload, tol)
                
                newCode = state2code(newState)
                A[newCode-1,code-1] += floatMotors*motor.Pattach / (1 + maxReattachLoc - minReattachLoc)
                D[newCode-1,code-1] += (new_x_c - shifted_x_c) * floatMotors * motor.Pattach / (1 + maxReattachLoc - minReattachLoc)
                
    diag = np.ravel(-1*A.sum(0)) #Column sum of sparse matrix A
    A.setdiag(diag) #Setting the diagonal to be sum of column without a for loop      
    
    # for ii in range(0,N):
    #     A[ii,ii] =  - np.sum(A[:,ii])          

    return A, D


#==============================================================================
# Set of functions to get the initial PDF of all configurations
'''
There are 3 functions in this section.
Take in number of motors, Load Force and motor parameters as inputs. 
More explanation to follow.
'''
#============================================================================== 
def initPDF(motor,m):
    numPos = int(2* np.floor(motor.l0/motor.ds) + 1)
    n = Nconfiguration_detach(numPos, m)
    
    stateVec = np.zeros((n,numPos + 2), dtype = 'float32')
    pdf = np.zeros(n,dtype = 'float32')
    
    for ii in range(1,n+1):
        state, Nm = code2state(ii, numPos)
        stateVec[ii-1,:numPos] = state
        
    for ii in range(0,n):
        ma = 0 #number of motors attached
        c = 0 #number of clusters
        s = 0 #spread
        mcvec = []
        
        for jj in range(0,numPos):
            sv = stateVec[ii][jj]
            if sv > 0: # i.e. motor is attached to MT at that location
                c += 1 #cluster found, update c by 1
                ma += sv #count the number of motors attached
                s = jj + 1 #spread = current location under consideration, i.e. jj
                mcvec.append(sv)
                
        nr = 1
        dr = 1
        
        for kk in range(0,c-1):
            nr = nr*np.math.factorial(ma)
            ma = ma - mcvec[c-1]
            dr = dr*np.math.factorial(int(mcvec[c-1]))
            
        mc = nr/dr
        deg = mc*(numPos - s + 1)
        stateVec[ii][numPos] = deg
    
    #Inserting the probabilities into the 2 extra columns  
    '''
    At least one motor assumed attached to the cargo, so the 0-motor attached 
    case is removed, here
    '''
    stateVec[0][numPos] = 0 
    normalizingSum = np.sum(stateVec,axis = 0)[-2] #Sum of 2nd last column
        
    for ii in range(0,n):
        stateVec[ii, numPos + 1] = (stateVec[ii,numPos]) / normalizingSum
        pdf[ii] = (stateVec[ii,numPos]) / normalizingSum

    return pdf, stateVec

def initPDFmod(motor,m,Fload):
    n = Nlocations(m, motor, Fload)
    
    '''
    Some doubt about what n should be here! look deeper into this.
    '''
    
#    n = int(2* np.floor(motor.l0/motor.ds) + 1)
    N = Nconfiguration_detach(n,m)
    pdf = np.zeros(N,dtype='float32')
    
    pdf_nz, stateVec = initPDF(motor,m)
    
    for ii in range(0,len(pdf_nz)):
        tempVar = np.zeros(n,dtype='int16') #Integer (-32768 to 32767) valued array
        tempVar[0:len(stateVec[ii][:-2])] = (stateVec[ii][:-2])
        code = state2code(tempVar)
        pdf[code-1] = pdf_nz[ii]
        
    return pdf

'''
next step, define initialize pdf. 
'''
def initializePDF(motor,m,Fload,Pa): 
    n = Nlocations(m, motor, Fload)
    numMotors = np.empty(0,dtype = 'int16')
    
    for nm in range(0,m+1):
        numMotors = np.append(numMotors, nm*np.ones(Nconfiguration(n,nm)))

    pdf = initPDFmod(motor,m,Fload)   
    pdf[0] = 1
    pdf0 = np.zeros(len(pdf))
    
    for nm in range(0,m+1):
        tempVar = pdf * (np.array([ii==nm for ii in numMotors]))
        tempVar = tempVar / sum(tempVar)
        value = (np.math.factorial(m)/(np.math.factorial(nm)*np.math.factorial(m-nm))*( (1-Pa)**(m-nm)*Pa**(nm) ))
        tempVar = tempVar * value
        pdf0 = pdf0 + tempVar
        
    pdf0 = pdf0 / sum(pdf0)
    
    return pdf0
    
#%% 
l0 = 110
Kel = 0.32e-3
ds = 8
Fs = 0.006
Pattach = 5
motor_type = 'kinesin'
Fload = 0.002
m = 3
ATP = 0.002
tol = 1e-6
motor = createMotor(l0, Kel, Fs, ATP, Pattach, motor_type)
A, D = getTransitionMatrix(motor, Fload, m, tol)
Aold = A.todense()
pdf = initializePDF(motor,m,Fload,1)

#%%









