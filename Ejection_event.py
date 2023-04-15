# -*- coding: utf-8 -*-
"""@author: bart_"""
######################################################################
import sys
import os
import numpy as np
import pandas as pd
import h5py as h5
import matplotlib.pyplot as plt

h = 0.6766
rhom = 8.634164473613977e-09
rhoc = 2.7753662724570803e-08

def radiic(rads,denz):
    a= np.nanmin(rads)
    b= np.nanmax(rads)
    global m
    rads = rads[~np.isnan(rads)]
    radius = np.sort(rads)
    while (b - a)/2.0 > 1E-12:
        midpoint = (a + b)/2.0
        partii = radius[radius <= midpoint]
        volume = (4.0/3.0) * np.pi * midpoint**3
        dens = m * partii.size / volume
        if dens < denz: # Increasing but below 0 case
            b = midpoint
        else:
            a = midpoint
    print(midpoint,m * partii.size)
    return midpoint,m * partii.size


np.random.seed(42069)
grow_radius = bool(sys.argv[3])

mass_fraction =    np.float64(sys.argv[1])     #%
ejection_radius =  27 #np.float64(sys.argv[2]) #[kpc]
files = np.array(os.listdir())
files = files[np.flatnonzero(np.core.defchararray.find(files,'output')!=-1)]
files = files[np.flatnonzero(np.core.defchararray.find(files,'.hdf5')!=-1)]
fil_num = np.int64(pd.Series(files).str[7:11].to_numpy()).max()

f = h5.File(f"output_{fil_num:04}.hdf5",'r')
coordinates = f["PartType1"]["Coordinates"]
#print(coordinates)
#coordinates = coordinates[~np.isnan(coordinates)].reshape(-1,3)
mass = f["PartType1"]["Masses"]
m = mass[0]
#mass = mass[~np.isnan(mass)]
velos = f["PartType1"]["Velocities"]
#velos = velos[~np.isnan(velos)].reshape(-1,3)
IDs = f["PartType1"]["ParticleIDs"]
#IDs = IDs[~np.isnan(IDs)]
total_mass = np.nansum(mass)
CoM = np.nanmean(coordinates,axis=0)
CoV = np.nanmean(velos,axis=0)
Radii = np.linalg.norm(coordinates-CoM,axis=1)
rvir,mvir = radiic(Radii,rhom*200)
#particles_total = np.sum(Radii<rvir)
particles_total = mvir / m
particles_to_eject = round(particles_total * mass_fraction/100 / 1.11)

if grow_radius:
	radii_sorted = np.argsort(Radii)
	mask = np.sort(radii_sorted[particles_to_eject:])

else:
	mask = Radii < ejection_radius
	particles_within = mask.sum()
	while particles_within < particles_to_eject:
		ejection_radius += 10
		mask = Radii < ejection_radius
		particles_within = mask.sum()
	particle_fraction = particles_to_eject/particles_within
	sample_mask = np.ones_like(Radii)
	eject_index = np.sort(np.random.choice(np.where(mask)[0], size=particles_to_eject, replace=False))
	sample_mask[eject_index] = False
	mask = np.where(sample_mask == 1)[0]

sample_c = coordinates[mask]
sample_v = velos[mask]
sample_m = mass[mask]
sample_IDs = IDs[mask]

#%% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#import the particle IC's
coords = sample_c
Mass = sample_m
PartID = sample_IDs
velos = sample_v

# dt = h5.vlen_dtype(np.dtype('int32'))
max_box = float(sys.argv[2])

#Make the .hdf5 file
with h5.File("IC_event.hdf5", "w") as fh5:
    head = fh5.create_group('Header')
    head.attrs.create('BoxSize', np.array([max_box, max_box, max_box]), dtype=np.float64)
    head.attrs.create('Dimension', np.array([3]), dtype=np.int32)
    head.attrs.create('Flag_Entropy_ICs', np.array([0, 0, 0, 0, 0, 0]), dtype=np.uint32)
    head.attrs.create('MassTable', np.array([0, 0, 0, 0, 0, 0]), dtype=np.float64)
    head.attrs.create('NumFilesPerSnapshot', 1, dtype=np.int32)
    head.attrs.create('NumPart_ThisFile', np.array([0,PartID.size,0,0,0,0]), dtype=np.uint32)
    head.attrs.create('NumPart_Total', np.array([0,PartID.size,0,0,0,0]), dtype=np.int32)
    head.attrs.create('NumPart_Total_HighWord', np.array([0,0,0,0,0,0]), dtype=np.int32)
    head.attrs.create('Redshift', np.array([0]), dtype=np.float64)
    head.attrs.create('Time', np.array([0]), dtype=np.float64)
    
    prt1 = fh5.create_group('PartType1')
    fh5.create_dataset('PartType1/Coordinates',data=coords,dtype=np.float64)
    fh5.create_dataset('PartType1/Masses',data=Mass,dtype=np.float16)
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

