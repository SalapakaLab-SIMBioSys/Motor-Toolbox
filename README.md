# Motor Toolbox
A toolbox for simulating intracellular cargo transport by a team of molecular motors.

## How to use
See motor_Example.py

To open the UserGuide, clone the Motor-Toolbox reppsitory and open the html files in any browser.

## Sample Code
```python
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
```

## Example Plots

### For Kinesin
![runLengthVsFload](https://user-images.githubusercontent.com/52796974/164310588-a65dae94-63aa-42cb-ba0d-18d0d267c2a7.png)
![velVsFload](https://user-images.githubusercontent.com/52796974/164310592-8d0b4a01-8bc1-41f2-b7dd-3def255be59c.png)


## Citation 
To cite SDA, please use the following publications:

- Materassi, Donatello, Subhrajit Roychowdhury, Thomas Hays, and Murti Salapaka. "An exact approach for studying cargo transport by an ensemble of molecular motors." BMC biophysics 6, no. 1 (2013): 1-18.


## Additional References

- Bhaban, Shreyas, Donatello Materassi, Mingang Li, Thomas Hays, and Murti Salapaka. "Interrogating emergent transport properties for molecular motor ensembles: A semi-analytical approach." PLoS computational biology 12, no. 11 (2016): e1005152.

- Shrivastava, Rachit, Sivaraj Sivaramakrishnan, and Murti V. Salapaka. "Cargo-motor interaction kinetics regulate myosin VI based transport." Biophysical Journal 121, no. 3 (2022): 402a.
