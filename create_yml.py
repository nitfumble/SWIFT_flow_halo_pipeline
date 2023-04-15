# -*- coding: utf-8 -*-
"""
@author: bart_
"""

import yaml
import sys
import numpy as np
import pandas as pd
np.random.seed(42069)
#check name yml
try:
    input_filename = str(sys.argv[1])
except (IndexError,TypeError,ValueError) as e:
    input_filename = 'isolated_galaxy.yml'
    
try:
    output_filename = 'IC_'+str(sys.argv[2])+'.hdf5'
except (IndexError,TypeError,ValueError) as e:
    output_filename = 'IC_Modelname.hdf5'

num_part = int(float(sys.argv[4]))

#Getting information on the dimensionallity:
max_box = float(sys.argv[3])
centre_box = max_box/2
input_filename2 = input_filename[:-4]+"_ejection"+input_filename[-4:]

#L_softening = (max_box/(num_part)**(1/3)/((num_part)**(1/3))**(1/3))/20
L_softening = ((max_box/2)/(num_part)**(1/3))/20
print(f'The Softening length is {L_softening} kpc \n')
Ejection_time = float(sys.argv[5]) *  3.15576000E16   #[s]

d = {    
# Define the system of units to use internally.
'InternalUnitSystem':{
  'UnitMass_in_cgs':                1.98848e+43     # 10^10 solar masses 
  ,'UnitLength_in_cgs':             3.08567758E21   # 1 kpc 
#  ,'UnitVelocity_in_cgs':          1E5             # km/s
  ,'UnitVelocity_in_cgs':           3.08567758E21   # kpc/s
  ,'UnitCurrent_in_cgs':            1               # Amperes
  ,'UnitTemp_in_cgs':               1               # Kelvin
  },

# Parameters for the self-gravity scheme
'Gravity':{
  'eta':                            0.025           # Constant dimensionless multiplier for time integration.
  ,'MAC':                           'geometric'
  ,'theta_cr':                      0.5             # Opening angle (Multipole acceptance criterion).
,'max_physical_DM_softening':     L_softening     # Physical softening length (in internal units).
  #,'comoving_DM_softening':     L_softening     # Physical softening length (in internal units).
#  ,'max_physical_baryon_softening': 0.01           # Physical softening length (in internal units).
},
# Parameters governing the time integration (Set dt_min and dt_max to the same value for a fixed time-step run.)
'TimeIntegration':{
  'time_begin':                     0.              # The starting time of the simulation (in internal units).
  ,'time_end':                      Ejection_time   #10     # The end time of the simulation (in internal units).
  ,'dt_min':                        3.15576000E10   #1e-6    # The minimal time-step size of the simulation (in internal units).
  ,'dt_max':                        3.15576000E14   #1e-2    # The maximal time-step size of the simulation (in internal units).
},
# Parameters governing the snapshots
'Snapshots':{
  'basename':                       'output'        # Common part of the name of output files
#  ,'time_first':                    0.             # (Optional) Time of the first output if non-cosmological time-integration (in internal units)
  ,'delta_time':                    3.15576000E15   # Time difference between consecutive outputs (in internal units)
},
#'Scheduler':{
#  'max_top_level_cells':            16
# }, 
# Parameters governing the conserved quantities statistics
'Statistics':{
  'delta_time':                     3.15576000E15   # Time between statistics output
#  ,'time_first':                    0.             # (Optional) Time of the first stats output if non-cosmological time-integration (in internal units)
},
# Parameters related to the initial conditions
'InitialConditions':{
   'file_name':                     output_filename     # The file to read
  ,'periodic':                      0                   # Are we running with periodic ICs?
}
}

# Parameters related to the particles in the given potential model
pos = {'Potential parameters':{
   'useabspos':                     0                   # Whether to use absolute position (1) or relative potential to centre of box (0)
  ,'position':                      [centre_box,centre_box,centre_box]       # Location of centre of potential with respect to centre of the box (internal units)
 # ,'epsilon':                       0.                 # No softening at the centre of the halo
 # ,'idealizeddisk':                 1                  # Run with an idealized galaxy disk
 # ,'h':                             0.7                # Reduced Hubble constant (value does not specify the used units!)
 # ,'vrot':                          200.               # Rotation speed of isothermal potential in internal units 
 # ,'scalelength':                   10.0               # Scale length of the potential
 # ,'diskfraction':                  0.0434370991372    # Disk mass fraction
 # ,'bulgefraction':                 0.00705852860979   # Bulge mass fraction
 # ,'timestep_mult':                 0.01               # Dimensionless pre-factor for the time-step condition, basically determines the fraction of the orbital time we use to do the time integration
 # ,'concentration':                 8.                 # Concentration of the Halo
 # ,'M_200':                         2.0e+12            # Virial mass (internal units)
}
}


with open(input_filename, 'w') as yaml_file:
    yaml.dump(d, yaml_file, default_flow_style=False)
    yaml.dump(pos, yaml_file, default_flow_style=None)
    # yaml.dump(pot, yaml_file, default_flow_style=False)
    
    
#%% From the example yml files

# #NFWPotential:
# NFW_MNPotential:
#   useabspos:        0          # 0 -> positions based on centre, 1 -> absolute positions 
#   position:         [0.,0.,0.] # Location of centre of isothermal potential with respect to centre of the box (if 0) otherwise absolute (if 1) (internal units)
#   timestep_mult:    0.01       # Dimensionless pre-factor for the time-step condition, basically determines the fraction of the orbital time we use to do the time integration
#   epsilon:          0.01       # Softening size (internal units)
#   concentration:    10.0       # concentration of the Halo
#   M_200:            150.0      # M200 of the galaxy disk
#   critical_density: 1.37E-8    # Critical density of the Universe in internal units
#   Mdisk:            3.0        # Disk mass (internal units)
#   Rdisk:            4.0        # Disk size (internal units)
#   Zdisk:            0.4704911  # Disk scale-height (internal units)


# # Hernquist potential parameters
# HernquistPotential:
#   useabspos:       0        # 0 -> positions based L_softeningon centre, 1 -> absolute positions 
#   position:        [0.,0.,0.]    # Location of centre of isothermal potential with respect to centre of the box (if 0) otherwise absolute (if 1) (internal units)
#   idealizeddisk:   1        # Run with an idealized galaxy disk
#   M200:            30.0   # M200 of the galaxy disk
#   h:               0.704    # reduced Hubble constant (value does not specify the used units!)
#   concentration:   7.1      # concentration of the Halo
#   diskfraction:              0.0434370991372   # Disk mass fraction
#   bulgefraction:              0.00705852860979  # Bulge mass fraction
#   timestep_mult:   0.01     # Dimensionless pre-factor for the time-step condition, basically determines the fraction of the orbital time we use to do the time integration
#   epsilon:         0.030      # Softening size (internal units)

last_file = int(d['TimeIntegration']['time_end']/d['Snapshots']['delta_time'])
d['Snapshots']['basename']          = 'eject'
#d['InitialConditions']['file_name'] = f'output_{last_file:04d}.hdf5'  
d['InitialConditions']['file_name'] = "IC_event.hdf5" 
#d['TimeIntegration']['time_begin']  = 3.15576000E17
d['TimeIntegration']['time_end']    = (4.418064E17-Ejection_time)

with open(input_filename2, 'w') as yaml_file:
    yaml.dump(d, yaml_file, default_flow_style=False)
    yaml.dump(pos, yaml_file, default_flow_style=None)
    # yaml.dump(pot, yaml_file, default_flow_style=False)
