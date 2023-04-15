# -*- coding: utf-8 -*-
"""@author: bart_"""
######################################################################
from numba import prange,config,jit,njit, float64
import matplotlib.pyplot as plt
from plotterfunc import *
import pandas as pd
import numpy as np
import h5py as h5
import imageio
import tables
import sys
import os

h = 0.6766
rho_m = 8.634164473613977e-09
rho_c = 2.7753662724570803e-08

config.THREADING_LAYER = 'omp'
cpuCount = os.cpu_count()
np.random.seed(42069)
print("Number of CPUs in the system:", cpuCount)
plot_bool = str(sys.argv[2])
if plot_bool == 'T':
	plot_bool = True
else:
	plot_bool = False


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


######################################################################
#%%  

paramfile = 'parameterfile.txt' 

d = {} #create dictionary
with open(paramfile) as f:
	for line in f:
		myline = line.split()
		if len(myline)>=2 and line[0]!='#': #make sure not blank line or comment
			(key,val) = myline[0:2] #extract first two items
			if key == 'Model':
				d[key] = val
			elif val == 'True':
				d[key] = True
			elif val == 'False':
				d[key] = False
			else:
				d[key] = float(val)

model = d['Model']
Grav_const = d['G']
params = d

data_modelname = str(sys.argv[1])
output_foldername = '../Figures/'
filename = 'haloplot' 
######################################################################
#%% 

if not os.path.exists('../Figures/images'):
	os.mkdir('../Figures/images')

files = os.listdir()
Nmax = np.flatnonzero(np.core.defchararray.find(files,'output')!=-1).size -3
Nmax_e =  np.flatnonzero(np.core.defchararray.find(files,'eject')!=-1).size -3
steps = 1
n_frames = 1
iterarray = np.arange(0, Nmax, steps)
iterarray_e = np.arange(0, Nmax_e, steps)

######################################################################
#%% INCL. EJECT


f = [h5.File("output_%04d.hdf5" % i, "r") for i in iterarray]
Numpartz = np.array([np.array(i["PartType1"]["Coordinates"].shape[0]) for i in f], dtype=np.float64)

coordinates_lists = np.array([np.concatenate((np.array(i["PartType1"]["Coordinates"]),np.array([np.nan]*3*(np.int64(d['n']-i["PartType1"]["Coordinates"].shape[0]))).reshape(-1,3)),axis=0) for i in f], dtype=np.float64)

masses = np.array([np.concatenate((np.array(i["PartType1"]["Masses"]),np.array([np.nan]*(np.int64(d['n']-i["PartType1"]["Masses"].shape[0]))).reshape(-1)),axis=0) for i in f], dtype=np.float64)

times = np.array([np.array(i["Header"].attrs["Time"][0], dtype=np.float) * 
          np.array(i["Units"].attrs["Unit time in cgs (U_t)"][0]) / 31557600.0e9 for i in f], dtype=np.float32)

velos_list = np.array([np.concatenate((np.array(i["PartType1"]["Velocities"]),np.array([np.nan]*3*(np.int64(d['n']-i["PartType1"]["Velocities"].shape[0]))).reshape(-1,3)),axis=0) for i in f])

#%% Parameters
combine = True 					#input('Combine the Density & Mass Profiles in 1 subplot? (Y/N)')
shift = np.array(f[0]["Header"].attrs["BoxSize"] / 2, dtype=np.float16) 	#float(input('What is the shift'))
#G = 4.299581e04 #??
figsize = 10

xmin = 10**.8
xmax = f[0]["Header"].attrs["BoxSize"][0]*10**.5
ymin0 = 10**-13.2
ymax0 = 10**-3.8
ymin1 = np.array(masses[0]).min()*10**-.2
ymax1 = np.array(masses[0]).sum()*10**.2
R_max = np.linalg.norm(np.array(f[0]["Header"].attrs["BoxSize"], dtype=np.float16))

[i.close() for i in f]
del f

for iter_part_e in np.array_split(iterarray_e,np.ceil(iterarray_e.size/500)):
	t_init = times[-1]
	f_e = [h5.File("eject_%04d.hdf5" % i, "r") for i in iter_part_e]

	Numpartz = np.concatenate((Numpartz,[np.array(i["PartType1"]["Coordinates"].shape[0]) for i in f_e]),axis=0)

	coordinates_lists = np.concatenate((coordinates_lists,[np.concatenate((np.array(i["PartType1"]["Coordinates"]),np.array([np.nan]*3*(np.int64(d['n']-i["PartType1"]["Coordinates"].shape[0]))).reshape(-1,3)),axis=0) for i in f_e]),axis=0)

	masses = np.concatenate((masses,[np.concatenate((np.array(i["PartType1"]["Masses"]),np.array([np.nan]*(np.int64(d['n']-i["PartType1"]["Masses"].shape[0]))).reshape(-1)),axis=0) for i in f_e]),axis=0)

	times = np.concatenate((times,t_init + [np.array(i["Header"].attrs["Time"][0]) * 
		      np.array(i["Units"].attrs["Unit time in cgs (U_t)"][0]) / 31557600.0e9 for i in f_e]),axis=0)

	velos_list = np.concatenate((velos_list,[np.concatenate((np.array(i["PartType1"]["Velocities"]),np.array([np.nan]*3*(np.int64(d['n']-i["PartType1"]["Velocities"].shape[0]))).reshape(-1,3)),axis=0) for i in f_e]),axis=0)

	[i.close() for i in f_e]
	del f_e

G = 4.51846772E-29
figsize = 10
min_d,max_d = 0,0
min_m,max_m = 0,0
R_min = 1
filenames1, filenames2 = [],[]

#Fill each list with nans up to n
CoM = np.nanmean(coordinates_lists,axis=1)
CoV = np.nanmean(velos_list,axis=1)
radius_particles = np.linalg.norm(np.array(coordinates_lists)-shift,axis=-1)
nbsamples = 100
Rmin = np.nanmin(radius_particles)
Rmax = np.nanmax(radius_particles) #np.linalg.norm(np.array(f[0]["Header"].attrs["BoxSize"]))
print(Rmax)
rsp = np.logspace(np.log10(Rmin), np.log10(Rmax), nbsamples)
max_radius = Rmax

rar = np.linspace(0,Rmax,1000)

# Methods to compute density profile
M = masses[0,0]
r = radius_particles
m = masses[0,0]
mass = m 

######################################################################
#%% CALCULATE DENSITY PROFILE

@njit(float64[:](float64[:], float64, float64[:], float64),fastmath=True)
def calc_radial_density_profile(R, m, br, Rmax):
	"""
	Calculates the radial density profile of particles with positions
	given by the array R, using a bin size dr and a maximum radius Rmax.
	Returns an array with the radial density profile values at each radial
	distance r_i from the origin, where r_i = (i + 0.5) * dr.
	"""
	# Compute the number of bins needed for the given parameters
	Nbins = br.size
	# Initialize an array to store the particle counts in each bin
	counts = np.zeros(Nbins, dtype=np.float64)
	# Loop over all particles and accumulate their counts in the appropriate bin
	for i in range(R.shape[0]):
		r = R[i]
		if r >= Rmax:
			continue
		counts[np.where(r > br)[0][-1]] += 1
	# Compute the volume of each shell
	shell_volumes = 4 * np.pi / 3 * (br[1:]**3 - br[:-1]**3)
	# Divide the counts by the shell volumes to get the radial density profile
	density_profile = m * counts[1:] / shell_volumes
	return density_profile
	
drad = 5
densities = np.empty(shape=((Nmax+Nmax_e+2,int(np.ceil(Rmax/drad))-1)))
binrange = np.logspace(-4,np.log10(Rmax),int(np.ceil(Rmax / drad)), dtype=np.float64)

for cc, radii in enumerate(radius_particles):
	densities[cc,:] = calc_radial_density_profile(radii, m, binrange, Rmax)

# remove empty bins to nan
densities = np.where(densities > 0, densities, np.nan)
densities[densities < 1E-14] = np.nan
densities[densities > 1E4] = np.nan
Dmin = np.nanmin(densities)
Dmax = np.nanmax(densities)

def find_radius_and_mass_at_density(R, m, denz):
	"""
	Uses the bisection method to find the radius at which a given density
	`rho_target` is achieved to within a tolerance `tol`, given an array of
	radii `R` and a maximum radius `Rmax`. Returns the radius and the mass
	enclosed within the radius.
	"""
	a= np.nanmin(R)
	b= np.nanmax(R)
	R = R[~np.isnan(R)]
	radius = np.sort(R)
	while (b - a)/2.0 > 1E-12:
		midpoint = (a + b)/2.0
		partii = radius[radius <= midpoint]
		volume = (4.0/3.0) * np.pi * midpoint**3
		dens = m * partii.size / volume
		if dens < denz: # Increasing but below 0 case
			b = midpoint
		else:
			a = midpoint

	#print(midpoint, m * partii.size)
	return midpoint, m * partii.size

#@njit((float64[:],float64,float64(float64[:,:],float64[:,:],float64[:,:], float64[:,:], float64), parallel=True, fastmath=True)
def bound_checker(m,v,c,r,v_com):
	v_net = (v-v_com)
	EK = 0.5*m*(v_net[:,0]**2+v_net[:,1]**2+v_net[:,2]**2)
	c = c[~np.isnan(c)].reshape(-1,3)
	potential = np.empty_like(c)
	for i in prange(c.shape[0]):
		potential[i] = np.sum(c-c[i],axis=0)/c.shape[0]
	EP = Grav_const*m**2*np.sqrt(potential[:,0]**2+potential[:,1]**2+potential[:,2]**2)
	bound_states = np.where(EK<EP,True,False)
	return np.array(bound_states),np.sum(EK),np.sum(EP)

######################################################################
#%% Calculating Overdensity Radii & Plotting
rm_200m_m_list = []
rm_200m_r_list = []
rm_500c_m_list = []
rm_500c_r_list = []
timestep_list = []
bound_perc_list = []
kin_list = []
pot_list = []
filenamez = np.array([])
Total_E = np.array([])
nframes = 1

f = h5.File("output_0000.hdf5", "r")

for index, radii in enumerate(radius_particles):
	framenr = index
	timestep_list.append(times[index])
	rm_500c = find_radius_and_mass_at_density(radii, m, 500*rho_c)
	rm_500c_r_list.append(rm_500c[0])
	rm_500c_m_list.append(rm_500c[1])
	rm_200m = find_radius_and_mass_at_density(radii, m, 200*rho_m)
	rm_200m_r_list.append(rm_200m[0])
	rm_200m_m_list.append(rm_200m[1])
	if plot_bool:
		nan_mask = np.array(~np.isnan(radii))	#Ignoring the nans!
		#print(f'Bound calculation: {index}',end='\r') 
		#bound,Kinetic_E,Potential_E = bound_checker(masses[index][:][nan_mask],velos_list[index][:][nan_mask],coordinates_lists[index][:][nan_mask]-CoM[index],radii[nan_mask],CoV[index])
		#kin_list.append(Kinetic_E)
		#pot_list.append(Potential_E)
		#Total_E = np.append(Total_E,np.array([Kinetic_E,Potential_E]),axis=0)
		bound = np.ones_like(coordinates_lists[index,:,0][nan_mask]) == 1
		Bound_perc = bound.sum()/bound.size*100
		bound_perc_list.append(Bound_perc)
		print('\t Bound %= '+"{:.2f}".format(Bound_perc),'% | @ t = '+str(np.round(times[index],1))+' Gyr',end='\r')
		filenamez = np.append(filenamez,np.array([output_foldername+'images/'+filename+str(framenr)]*nframes))
		c = densities[index] > 0
		dens = np.compress(c, densities[index])
		rs_i = np.compress(c, np.diff(binrange))


######################################################################
#%% CALCULATE THE 200M AND 500C MASS & RADII
#dens_c_now = 1.32934445e-8
#Might be usefull:
#from colossus.halo import mass_so
#from hmf import *

# Crit dens univ now: 9.47 x 10^-27 kg m^-3
# Crit dens univ now: 132.934445 solar mass (kpc^(-3))
# Mean dens univ now: 9.9×10-27 kg m-3
# Mean dens univ now: 146.22789 solar mass (kpc^(-3))


######################################################################
#%% PLOTTING THE FIGURE(S)

		plot3d_2scat_prof(
	#Bound:
		x=coordinates_lists[index,:,0][nan_mask][bound],
		y=coordinates_lists[index,:,1][nan_mask][bound],
		z=coordinates_lists[index,:,2][nan_mask][bound],
	#Unbound:
		x1=coordinates_lists[index,:,0][nan_mask][bound != True],
		y1=coordinates_lists[index,:,1][nan_mask][bound != True],
		z1=coordinates_lists[index,:,2][nan_mask][bound != True],
	#Sphereical Overdensity:
		centre=CoM[index,:],
		#(shift*2,shift*2),
		radi_200m=rm_200m[0],
		m200m=rm_200m[1],
		radi_500c=rm_500c[0],
		m500c=rm_500c[1],
	#Plotting Info:
		lims=[-0.1*np.array(f["Header"].attrs["BoxSize"])[0], 1.1*np.array(f["Header"].attrs["BoxSize"])[0]],
		folder=output_foldername+'images/',
		name=filename+str(framenr),
		title='Spatial distribution of halo particles',
		time=times[index],
		fbound=Bound_perc,
		rs=rs_i,
		dens=dens,
		Rlims=(0.2*Rmin,Rmax),
		Denslims=(Dmin,Dmax),
		dens_model=[])

#Bounds = np.array(bound_perc_list)
timestep = np.array(timestep_list).ravel()

plot_evo(timestep,np.array(rm_500c_r_list),r'R$_{500c}$',np.array(rm_200m_r_list),r'R$_{200m}$',f'Evolution of the overdensity radii ({data_modelname})','Time [Gyr]',r'Radius $h^{-1}$ [kpc]',output_foldername)
plot_evo(timestep,np.array(rm_500c_m_list),r'M$_{500c}$',np.array(rm_200m_m_list),r'M$_{200m}$',f'Evolution of the overdensity masses ({data_modelname})','Time [Gyr]',r'Mass $\times$ 10$^{10}$ $h^{-1}$ [M$_{\odot}$]',output_foldername)

plot_evo(timestep,np.array(rm_500c_r_list)/np.array(rm_500c_r_list)[0],r'R$_{500c}$',np.array(rm_200m_r_list)/np.array(rm_200m_r_list)[0],r'R$_{200m}$',f'Ratio Evolution of the overdensity radii ({data_modelname})','Time [Gyr]','Ratio',output_foldername)
plot_evo(timestep,np.array(rm_500c_m_list)/np.array(rm_500c_m_list)[0],r'M$_{500c}$',np.array(rm_200m_m_list)/np.array(rm_200m_m_list)[0],r'M$_{200m}$',f'Ratio Evolution of the overdensity masses ({data_modelname})','Time [Gyr]','Ratio',output_foldername)

#sys.exit()

######################################################################
#% Build GIF        
if plot_bool:
	print('creating gif density profile\n')
	with imageio.get_writer(output_foldername+f'Overdensity region plot ({data_modelname}).gif', mode='I') as writer:
		for fn in filenamez:
			print('Loading: '+str(fn),end='\r')
			ims = imageio.imread(str(fn)+'.png')
			writer.append_data(ims)
			os.remove(str(fn)+'.png')

	print('gif complete\n')
	print('Done!')


dataset = pd.DataFrame({
		'Times': np.array(timestep_list), 
		'radi_200m': np.array(rm_200m_r_list), 
		'm200m': np.array(rm_200m_m_list), 
		'radi_500c': np.array(rm_500c_r_list), 
		'm500c': np.array(rm_500c_m_list), 
		#'Total_Energy': np.sum(Total_E,axis=1),
		#'Potential_Energy': np.array(Total_E[:,1]), 
		#'Kinetic_Energy': np.array(Total_E[:,0]), 
		'Number_of_particles': Numpartz,
		#'Bound_perc': Bounds,
		})
dataset.to_csv(f'../Results/Halo_vision ({data_modelname}).csv',index=False)


############################################################################################

