# -*- coding: utf-8 -*-
"""@author: bart"""

import sys
import os
import numpy as np
import pandas as pd
import h5py as h5

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
    input_mass = 0.00014242566 #10^10 M_sun

#import the particle IC's
df = pd.read_csv(input_filename)
ID =  df[df.columns[0]]
x =  df.x.values
y =  df.y.values
z =  df.z.values
vx =  df.vx.values
vy =  df.vy.values
vz =  df.vz.values

coords = np.stack((x,y,z),axis=1)
Mass = np.array([input_mass]*x.shape[0])
PartID = ID.values + 1
velos = np.stack((vx,vy,vz),axis=1)

# dt = h5.vlen_dtype(np.dtype('int32'))
max_box = float(sys.argv[3])

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
    fh5.create_dataset('PartType1/Coordinates',data=coords,dtype=np.float32)
    fh5.create_dataset('PartType1/Masses',data=Mass,dtype=np.float16)
    fh5.create_dataset('PartType1/ParticleIDs',data=PartID,dtype=np.int16)
    fh5.create_dataset('PartType1/Velocities',data=velos,dtype=np.float32)
    
    # prt4 = fh5.create_group('PartType4')
    # fh5.create_dataset('PartType4/Coordinates',track_order=True)
    # fh5.create_dataset('PartType4/Masses',track_order=True)
    # fh5.create_dataset('PartType4/ParticleIDs',track_order=True)
    # fh5.create_dataset('PartType4/Velocities',track_order=True)
    
    
    runt = fh5.create_group('RuntimePars')
    runt.attrs.create('PeriodicBoundariesOn', 0, dtype=np.int32)
    
    unts = fh5.create_group('Units')
    unts.attrs.create('Unit current in cgs (U_I)', np.array([1]), dtype=np.float64)
    unts.attrs.create('Unit length in cgs (U_L)', np.array([3.08567758E+21]), dtype=np.float64)
    unts.attrs.create('Unit mass in cgs (U_M)', np.array([1.98848e+43]), dtype=np.float64)
    unts.attrs.create('Unit temperature in cgs (U_T)', np.array([1]), dtype=np.float64)
    unts.attrs.create('Unit time in cgs (U_t)', np.array([1]), dtype=np.float64)
    #unts.attrs.create('UnitVelocity_in_cgs', np.array([3.08567758E+21]), dtype=np.float64)
# h5.File.close(fh5)
#'UnitVelocity_in_cgs':           3.08567758E21   # kpc/s

#%% NOT THE CODE
#Just to check with other known files
# a = h5.File("3e11-star-only-DM-halo-galaxy(2).hdf5")
# list(a)

# alist = [i for i in a] 
# blist = [a[i] for i in alist]
# clist = [i for i in blist]

# 'Header'
# []

# 'PartType1'
# ['Coordinates', 'Masses', 'ParticleIDs', 'Velocities']

# 'PartType4'
# ['Coordinates', 'Masses', 'ParticleIDs', 'Velocities']

# 'RuntimePars'
# []

# 'Units'
# []



# b=a['PartType1']
# list(b)
# list(b['Coordinates'])[:5]
# list(b['Masses'])[:5]
# list(b['ParticleIDs'])[:5]
# list(b['Velocities'])[:5]




