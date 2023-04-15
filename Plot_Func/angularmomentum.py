#!/usr/bin/env python
###############################################################################
# This file is part of SWIFT.
# Copyright (c) 2018 Folkert Nobels (nobels@strw.leidenuniv.nl)
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published
# by the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
################################################################################
#
# Run with: python angularmomentum.py Figures/ NFW
#
################################################################################
import numpy as np
import pandas as pd
import h5py
import matplotlib.pyplot as plt
import scipy.optimize as sco
import sys
import os 
np.random.seed(42069)

#PLOT PARAMS

#  Set the font size for axis labels and tick labels
plt.rcParams['font.size'] = 18

# Set the font family for all text in the plot
plt.rcParams['font.family'] = 'serif'

# Set the figure size to 6 x 4 inches
plt.rcParams['figure.figsize'] = [10, 8]

# Set the linewidth for lines in the plot
plt.rcParams['lines.linewidth'] = 1.5

# Set the color cycle for multiple lines in the same plot
plt.rcParams['axes.prop_cycle'] = plt.cycler('color', ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf', '#a55194'])

# Set the tick direction to 'in'
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'

# Set the tick length to 4 points
plt.rcParams['xtick.major.size'] = 4
plt.rcParams['ytick.major.size'] = 4

# Set the number of minor ticks between major ticks to 5
plt.rcParams['xtick.minor.visible'] = True
plt.rcParams['ytick.minor.visible'] = True
plt.rcParams['xtick.minor.size'] = 2
plt.rcParams['ytick.minor.size'] = 2
plt.rcParams['xtick.minor.width'] = 0.5
plt.rcParams['ytick.minor.width'] = 0.5
plt.rcParams['xtick.minor.pad'] = 2.0
plt.rcParams['ytick.minor.pad'] = 2.0
plt.rcParams['xtick.minor.top'] = True
plt.rcParams['ytick.minor.right'] = True
plt.rcParams['xtick.minor.bottom'] = True
plt.rcParams['ytick.minor.left'] = True

# Set the default dpi of figures to 150
plt.rcParams['figure.dpi'] = 150

# Set the default file format to PDF
plt.rcParams['savefig.format'] = 'pdf'


files = os.listdir()
Nmax =  np.flatnonzero(np.core.defchararray.find(files,'output')!=-1).size -3
steps = 1
angmomcomp = True

iterarray = np.arange(0, Nmax, steps)
Lxtot = np.zeros(len(iterarray))
Lytot = np.zeros(len(iterarray))
Lztot = np.zeros(len(iterarray))
Ltot = np.zeros(len(iterarray))
time_array = np.zeros(len(iterarray))

output_foldername = str(sys.argv[1])
model = str(sys.argv[2])
data_foldername = 'SwiftRUN_'+model+'/'
#data_foldername = str(sys.argv[3])


for i in iterarray:
#    f = h5py.File(data_foldername+"output_%04d.hdf5" % i, "r")
    f = h5py.File("output_%04d.hdf5" % i, "r")
    
    boxsize = f["Header"].attrs["BoxSize"] / 2

    time_array[int(i / steps)] = f["Header"].attrs["Time"]

    particles = f["PartType1"]
    coordinates = particles["Coordinates"][:, :]
    velocities = particles["Velocities"][:, :]
    masses = particles["Masses"][:]

    R = (
        (coordinates[:, 0] - boxsize[0]) ** 2 + (coordinates[:, 1] - boxsize[1]) ** 2
    ) ** 0.5
    X = np.abs(coordinates[:, 0] - boxsize[0])
    Y = np.abs(coordinates[:, 1] - boxsize[1])
    Z = np.abs(coordinates[:, 2] - boxsize[2])

    vx = velocities[:, 0]
    vy = velocities[:, 1]
    vz = velocities[:, 2]

    Lx = (Y * vz - Z * vy) * masses
    Ly = (Z * vx - X * vz) * masses
    Lz = (X * vy - Y * vx) * masses

    L = (Lx ** 2 + Ly ** 2 + Lz ** 2) ** 0.5

    Lxtot[int(i / steps)] = np.sum(Lx)
    Lytot[int(i / steps)] = np.sum(Ly)
    Lztot[int(i / steps)] = np.sum(Lz)
    Ltot[int(i / steps)] = np.sum(L)
    #print(Ltot[int(i / steps)])

Nmax_e =  np.flatnonzero(np.core.defchararray.find(files,'eject')!=-1).size -3
steps_e = 1
angmomcomp_e = True

iterarray_e = np.arange(0, Nmax_e, steps_e)
Lxtot_e = np.zeros(len(iterarray_e))
Lytot_e = np.zeros(len(iterarray_e))
Lztot_e = np.zeros(len(iterarray_e))
Ltot_e = np.zeros(len(iterarray_e))
time_array_e = np.zeros(len(iterarray_e))


for i in iterarray_e:
#    f = h5py.File(data_foldername+"eject_%04d.hdf5" % i, "r")
    f = h5py.File("eject_%04d.hdf5" % i, "r")
    boxsize = f["Header"].attrs["BoxSize"] / 2
    time_array_e[int(i / steps)] = time_array[-1]+ f["Header"].attrs["Time"]

    particles = f["PartType1"]
    coordinates = particles["Coordinates"][:, :]
    velocities = particles["Velocities"][:, :]
    masses = particles["Masses"][:]

    R = (
        (coordinates[:, 0] - boxsize[0]) ** 2 + (coordinates[:, 1] - boxsize[1]) ** 2
    ) ** 0.5
    X = np.abs(coordinates[:, 0] - boxsize[0])
    Y = np.abs(coordinates[:, 1] - boxsize[1])
    Z = np.abs(coordinates[:, 2] - boxsize[2])

    vx = velocities[:, 0]
    vy = velocities[:, 1]
    vz = velocities[:, 2]

    Lx = (Y * vz - Z * vy) * masses
    Ly = (Z * vx - X * vz) * masses
    Lz = (X * vy - Y * vx) * masses

    L = (Lx ** 2 + Ly ** 2 + Lz ** 2) ** 0.5

    Lxtot_e[int(i / steps)] = np.sum(Lx)
    Lytot_e[int(i / steps)] = np.sum(Ly)
    Lztot_e[int(i / steps)] = np.sum(Lz)
    Ltot_e[int(i / steps)] = np.sum(L)

iterarray = np.concatenate((iterarray,iterarray_e))
Lxtot = np.concatenate((Lxtot,Lxtot_e))
Lytot = np.concatenate((Lytot,Lytot_e))
Lztot = np.concatenate((Lztot,Lztot_e))
Ltot = np.concatenate((Ltot,Ltot_e))
time_array = np.concatenate((time_array,time_array_e))/3.15576000E16

fig0 = plt.figure()
if angmomcomp:
    plt.plot(time_array, Lxtot, label="Lx total")
    plt.plot(time_array, Lytot, label="Ly total")
    plt.plot(time_array, Lztot, label="Lz total")
#print(time_array)    
plt.plot(time_array, Ltot, label="L total")
plt.title('Angular Momentum (Model: '+model+')')
plt.xlabel("Time [Gyr]")
plt.ylabel(r"Angular momentum [$10^{10}$ M$_{\odot}$ kpc$^{2}$ s$^{-1}$]")
plt.legend()
plt.savefig(output_foldername+"Angular_momentum (Model "+model+").png")
#plt.show()
#plt.pause(0.0001)
plt.close(fig0)

fig1 = plt.figure()
#time_array[-1] = 2.0
if angmomcomp:
    plt.plot(time_array, Lxtot / Lxtot[0] - 1, label="Lx total")
    plt.plot(time_array, Lytot / Lytot[0] - 1, label="Ly total")
    plt.plot(time_array, Lztot / Lztot[0] - 1, label="Lz total")
#print(time_array)    
plt.plot(time_array, Ltot / Ltot[0] - 1, label="L total")
plt.title('Ratio Current and Zero Angular Momentum (Model: '+model+')')
plt.xlabel("Time [Gyr]")
plt.ylabel("Ratio")
plt.legend()
plt.savefig(output_foldername+"Ratio_angular_momentum (Model "+model+").png")
#plt.show()
#plt.pause(0.0001)
plt.close(fig1)

fig2 = plt.figure()
plt.title('Fractional Change of Total Angular Momentum (Model: '+model+')')
plt.semilogy(time_array, np.absolute(Ltot / Ltot[0] - 1))
plt.xlabel("Time [Gyr]")
plt.ylabel("Fraction")
plt.savefig(output_foldername+"Fractional_change_angular_momentum (Model "+model+").png")
#plt.show()
#plt.pause(0.0001)
plt.close(fig2)

dataset = pd.DataFrame({'Time': time_array, 'Total_angmom': Ltot, 'Total_x_angmom': Lxtot, 'Total_y_angmom': Lytot, 'Total_z_angmom': Lztot})
dataset.to_csv(f'../Results/Angular_momentum ({model}).csv',index=False)
