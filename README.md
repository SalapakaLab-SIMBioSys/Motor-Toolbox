# Motor Toolbox
A toolbox for simulating intracellular cargo transport by a team of molecular motors.

## How to use
See motor_Example.py for details.

To open the UserGuide, clone the Motor-Toolbox reppsitory and open the html files in any browser.

## Sample Code
```python
#==============================================================================
# Importing packages and motor toolbox
#==============================================================================
import numpy as np
from motor import *

#==============================================================================
# Defining motor, environment, and simulation parameters
#==============================================================================

l0 = 110 #rest length of motor in nm
Kel = 0.32e-3 #elasticity coefficient of motor in nN/nm
ds = 8 #step size of motor in nm
Fs = 0.006 #stall force of motor in nN
Pattach = 5 #reattachment rate of motor in per second
motor_type = 'kinesin'#motor type, 'kinsin', 'myosin', 'dynein'
Fload = 0.002 #load force on cargo in nN
tol = 1e-6 #tolerance for calculating equilibrium position of cargo
ATP = 0.002 #ATP concentration in M
K_Pd = 0.04 #constant with which stepping rate is multiplied to get detachment rate
Pback = 2 #detachment rate of motor after stall in per second
motor = createMotor(l0, Kel, Fs, ATP, Pattach, motor_type) #object of class motorClass 
mrange = range(2,4) #range of motors in an ensemble user wants to simulate
Frange = np.arange(0.002,0.012,0.0005) #load force range
Tend = 10 #total time horizon for a simulation in seconds

#==============================================================================
# Getting runlength and velocity plots
#==============================================================================

fig1 = getVelVsFload(motor, mrange, Frange, save = False, tol = 1e-6)
fig2 = getRunVsFload(motor, mrange, Frange, Tend, save = False, tol = 1e-6)
```

## Example Plots

### For Kinesin
<img src="https://user-images.githubusercontent.com/52796974/164310588-a65dae94-63aa-42cb-ba0d-18d0d267c2a7.png" width="400"> !<img src="https://user-images.githubusercontent.com/52796974/164310592-8d0b4a01-8bc1-41f2-b7dd-3def255be59c.png" width="400">


## Citation 
To cite Motor Toolbox, please use the following publications:

- Materassi, Donatello, Subhrajit Roychowdhury, Thomas Hays, and Murti Salapaka. "An exact approach for studying cargo transport by an ensemble of molecular motors." BMC biophysics 6, no. 1 (2013): 1-18.


## Additional References

- Bhaban, Shreyas, Donatello Materassi, Mingang Li, Thomas Hays, and Murti Salapaka. "Interrogating emergent transport properties for molecular motor ensembles: A semi-analytical approach." PLoS computational biology 12, no. 11 (2016): e1005152.

- Shrivastava, Rachit, Sivaraj Sivaramakrishnan, and Murti V. Salapaka. "Cargo-motor interaction kinetics regulate myosin VI based transport." Biophysical Journal 121, no. 3 (2022): 402a.
