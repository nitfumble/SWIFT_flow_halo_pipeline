# -*- coding: utf-8 -*-
"""@author: bart"""
import pandas as pd
import numpy as np
import h5py as h5
import yaml
import sys
import os

np.random.seed(42069)
try:
    os.remove("IC_"+str(sys.argv[2])+".hdf5")
except (FileNotFoundError,TypeError,ValueError) as e:
    pass

try:
    input_filename = str(sys.argv[1])
except (IndexError,TypeError,ValueError) as e:
    input_filename = 'outputfile.csv'

try:
   input_model = str(sys.argv[2])
except (IndexError,TypeError,ValueError) as e:
    input_model = 'Model'

try:
   input_mass = float(sys.argv[4])
except (IndexError,TypeError,ValueError) as e:
    input_mass = 0.0115 #[10^10 M_sun]  = M_MW

#import the particle IC's
df = pd.read_csv(input_filename)
ID =  df[df.columns[0]]
x =  df.x.to_numpy(dtype=np.float64)
y =  df.y.to_numpy(dtype=np.float64)
z =  df.z.to_numpy(dtype=np.float64)
vx =  df.vx.to_numpy(dtype=np.float64)
vy =  df.vy.to_numpy(dtype=np.float64)
vz =  df.vz.to_numpy(dtype=np.float64)

coords = np.stack((x,y,z),axis=1)
Mass = np.array([input_mass]*x.shape[0],dtype=np.float64)
PartID = ID.values + 1
velos = np.stack((vx,vy,vz),axis=1)

# dt = h5.vlen_dtype(np.dtype('int32'))
max_box = float(sys.argv[3])
num_part = int(float(sys.argv[5]))
centre_box = max_box/2

#Make the .hdf5 file
with h5.File("IC_"+input_model+".hdf5", "w") as fh5:

    head = fh5.create_group('Header')
    head.attrs.create('BoxSize', np.array([max_box, max_box, max_box]), dtype=np.float64)
    head.attrs.create('Dimension', np.array([3]), dtype=np.int32)
    head.attrs.create('Flag_Entropy_ICs', np.array([0, 0, 0, 0, 0, 0]), dtype=np.uint32)
    head.attrs.create('MassTable', np.array([0, 0, 0, 0, 0, 0]), dtype=np.float64)
    head.attrs.create('NumFilesPerSnapshot', 1, dtype=np.int32)
    head.attrs.create('NumPart_ThisFile', np.array([0,ID.size,0,0,0,0]), dtype=np.uint32)
    head.attrs.create('NumPart_Total', np.array([0,ID.size,0,0,0,0]), dtype=np.int32)
    head.attrs.create('NumPart_Total_HighWord', np.array([0,0,0,0,0,0]), dtype=np.int32)
    head.attrs.create('Redshift', np.array([0]), dtype=np.float64)
    head.attrs.create('Time', np.array([0]), dtype=np.float64)
    
    prt1 = fh5.create_group('PartType1')
    fh5.create_dataset('PartType1/Coordinates',data=coords,dtype=np.float64)
    fh5.create_dataset('PartType1/Masses',data=Mass,dtype=np.float64)
    fh5.create_dataset('PartType1/ParticleIDs',data=PartID,dtype=np.int16)
    fh5.create_dataset('PartType1/Velocities',data=velos,dtype=np.float64)
    
    runt = fh5.create_group('RuntimePars')
    runt.attrs.create('PeriodicBoundariesOn', 0, dtype=np.int16)
    
    unts = fh5.create_group('Units')
    unts.attrs.create('Unit current in cgs (U_I)', np.array([1]), dtype=np.float64)
    unts.attrs.create('Unit length in cgs (U_L)', np.array([3.08567758E+21]), dtype=np.float64)
    unts.attrs.create('Unit mass in cgs (U_M)', np.array([1.98848e+43]), dtype=np.float64)
    unts.attrs.create('Unit temperature in cgs (U_T)', np.array([1]), dtype=np.float64)
    unts.attrs.create('Unit time in cgs (U_t)', np.array([1]), dtype=np.float64)

L_softening = (max_box/(num_part)**(1/3)/((num_part)**(1/3))**(1/3))/20
#L_softening = (max_box/(num_part)**(1/3))/20

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
},
# Parameters governing the time integration (Set dt_min and dt_max to the same value for a fixed time-step run.)
'TimeIntegration':{
  'time_begin':                     0.              # The starting time of the simulation (in internal units).
  ,'time_end':                      1.893456E18   #10     # The end time of the simulation (in internal units).
  ,'dt_min':                        3.15576000E10   #1e-6    # The minimal time-step size of the simulation (in internal units).
  ,'dt_max':                        3.15576000E14   #1e-2    # The maximal time-step size of the simulation (in internal units).
},
# Parameters governing the snapshots
'Snapshots':{
  'basename':                       'IC'        # Common part of the name of output files
  ,'delta_time':                    1.893456E18   # Time difference between consecutive outputs (in internal units)
},
'Statistics':{
  'delta_time':                     4.73364E16   # Time between statistics output
},
# Parameters related to the initial conditions
'InitialConditions':{
   'file_name':                     "IC_"+input_model+".hdf5"     # The file to read
  ,'periodic':                      0                   # Are we running with periodic ICs?
}
}

# Parameters related to the particles in the given potential model
pos = {'Potential parameters':{
   'useabspos':                     0                   # Whether to use absolute position (1) or relative potential to centre of box (0)
  ,'position':                      [centre_box,centre_box,centre_box]       # Location of centre of potential with 
}
}


with open("Relax.yml", 'w') as yaml_file:
    yaml.dump(d, yaml_file, default_flow_style=False)
    yaml.dump(pos, yaml_file, default_flow_style=None)
