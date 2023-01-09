# -*- coding: utf-8 -*-
"""@author: bart_"""
######################################################################
from colossus.cosmology import cosmology
from scipy.optimize import curve_fit
from numba.np.ufunc import parallel
import matplotlib.pyplot as plt
from plotterfunc import *
from numba import prange,config,njit,vectorize,set_num_threads,get_num_threads
import numpy as np
import h5py as h5
import imageio
import sys
import os
config.THREADING_LAYER = 'tbb'
cpuCount = os.cpu_count()

print("Number of CPUs in the system:", cpuCount)
#set_num_threads(4)
print("Number of CPUs used:", get_num_threads())
######################################################################
#%% 

cosmo = cosmology.setCosmology('planck18')

def model_NFW(r,M_tot):
# Read the params from parameterfile.txt!
	rs=params['r_s']
	
	
	if params['truncate']:
		rc=params['r_cut']
	return
	
#def model_NFWX(): ICs not availible yet!
#	return

def model_Hernquist():
	return

def model_King():
	return

def model_Einasto():
	return
	
	
######################################################################
#%% 

paramfile = 'SwiftRUN_'+str(sys.argv[2])+'/parameterfile.txt' 

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
m = d['m']
Grav_const = d['G']
params = d

data_foldername = 'SwiftRUN_'+model+'/'
output_foldername = str(sys.argv[1])
filename = 'haloplot' 
######################################################################
#%% 

if not os.path.exists('Figures/images'):
	os.mkdir('Figures/images')

Nmax = 100
steps = 1
n_frames = 1

endp = 100# np.int8(input('How many files?'))
strings = [str(i) for i in np.arange(0,endp+1)]

nams = np.char.zfill(strings,3)
#fnames = [f'output_0{i}.hdf5' for i in nams]
#f = [h5.File(f, "r") for f in fnames]

iterarray = np.arange(0, Nmax + 1, steps)
f = [h5.File(data_foldername+"output_%04d.hdf5" % i, "r") for i in iterarray]

coordinates_lists = [i["PartType1"]["Coordinates"] for i in f]
masses = [i["PartType1"]["Masses"] for i in f]
times = [i["Header"].attrs["Time"][0] * 
          i["Units"].attrs["Unit time in cgs (U_t)"][0] / 31557600.0e9 for i in f]

velos_list = [i["PartType1"]["Velocities"] for i in f]

#%% Parameters
combine = True 					#input('Combine the Density & Mass Profiles in 1 subplot? (Y/N)')
shift = np.array(f[0]["Header"].attrs["BoxSize"] / 2) 	#float(input('What is the shift'))
#G = 4.299581e04 #??
figsize = 10

xmin = 10**.8
xmax = f[0]["Header"].attrs["BoxSize"][0]*10**.5
ymin0 = 10**-13.2
ymax0 = 10**-3.8
ymin1 = np.array(masses[0]).min()*10**-.2
ymax1 = np.array(masses[0]).sum()*10**.2


#if combine == 'Y':
#    combine = True
#else:
#    combine = False

G = 4.299581e04
figsize = 10
min_d,max_d = 0,0
min_m,max_m = 0,0
R_min = 1
R_max = np.linalg.norm(np.array(f[0]["Header"].attrs["BoxSize"]))
filenames1, filenames2 = [],[]
#Fill each list with nans up to n
filler = [np.int16(d['n']-np.array(fill).shape[0]) for fill in coordinates_lists]
new_coords_list = np.array([np.concatenate((np.array(C),np.array([np.nan]*3*filler[x]).reshape(filler[x],3)),axis=0) for x,C in enumerate(coordinates_lists)])
CoM = np.nanmean(new_coords_list,axis=1)
new_velos_list = np.array([np.concatenate((np.array(C),np.array([np.nan]*3*filler[x]).reshape(filler[x],3)),axis=0) for x,C in enumerate(velos_list)])
new_masses = np.array([np.concatenate((np.array(C),np.array([np.nan]*filler[x]).reshape(filler[x])),axis=0) for x,C in enumerate(masses)])
CoV = np.nanmean(new_velos_list,axis=1)
radius_particles = np.linalg.norm(np.array(new_coords_list)-shift,axis=-1)
#radius_particles = np.array([np.linalg.norm(np.array(C)-shift,axis=-1) for C in new_coords_list])
#max_radius = radius_particles.max(axis=0)

CoV_inner_particles = np.nanmean(np.where(radius_particles[:,:,np.newaxis]<40,new_velos_list,np.nan),axis=1)

nbsamples = 100
Rmin = np.nanmin(radius_particles)
Rmax = np.nanmax(radius_particles) #np.linalg.norm(np.array(f[0]["Header"].attrs["BoxSize"]))
rsp = np.logspace(np.log10(Rmin), np.log10(Rmax), nbsamples)
max_radius = Rmax
sort_r_particles = np.sort(radius_particles,axis=-1)
dens_arr = m/(4/3*np.pi*sort_r_particles**3)

	
rar = np.linspace(0,Rmax,1000)

#check = np.array([(sort_r_particles[sort_r_particles<r].size*m)/(4/3*np.pi*r**3) for r in rar])

rho_c = 139.876577E-10 	#10^10 solar mass (kpc^(-3))
Ω_m = 0.25
rho_m = Ω_m*rho_c	#10^10 solar mass (kpc^(-3))

# Methods to compute density profile
M = new_masses[0,0]
r = radius_particles
#mass = new_masses
mass = m 

"""#mass_shells = np.array([((r < R) * mass).sum(axis=1) for R in rar])
#mass_diff = np.diff(mass_shells,axis=0)
#R_diff = np.diff(rar)

q = np.array(mass_diff / R_diff[:, np.newaxis] / (4.0 * np.pi * rar[1:, np.newaxis] ** 2))

q_mask = np.where(q > 0,q,np.nan)
i0c, i1c = np.where(q_mask < 500*rho_c)
i0m, i1m = np.where(q_mask < 200*rho_m)

p1c = np.unique(i1c,return_index=True)
p0c = i0c[p1c[1]]
p1c = p1c[0]
posc = np.array([p0c,p1c])
r_500c= rar[p0c]

p1m = np.unique(i1m,return_index=True)
p0m = i0m[p1m[1]]
p1m = p1m[0]
posm = np.array([p0m,p1m])
r_200m= rar[p0m]"""

######################################################################
#%% CALCULATE DENSITY PROFILE


def mass_ins(R):
	return ((radii < R) * mass).sum(axis=-1)
	
def density(R):
	return np.diff(mass_ins_vect(R)) / np.diff(R) / (4.0 * np.pi * R[1:] ** 2)

# Methods to compute density profile

mass_ins_vect = np.vectorize(mass_ins)
rs = rsp[1:]
densities = np.empty(shape=((Nmax+1,rs.shape[0])))
for cc, radii in enumerate(radius_particles):
	densities[cc,:] = density(rsp)

# remove empty bins to nan
densities = np.where(densities > 0, densities, np.nan)
Dmin = np.nanmin(densities)
Dmax = np.nanmax(densities)

#print(r_500c,r_200m)
#print(r.min())
#print(np.sum(r<r_500c[:,np.newaxis]*m,axis=1))
#print(np.sum(r<r_200m[:,np.newaxis]*m,axis=1))
#plt.loglog(rar[1:,np.newaxis],q)
#plt.show()

def dense(R2,Rlim):
	mass = m*np.sum(R2<Rlim)
	#print(mass/m,R2.size)
	dens = mass/(4/3*np.pi*Rlim**3)
	return dens,mass

def bisec(a, b, R2, val, tol):
	if not (dense(R2,a)[0] > val and dense(R2,b)[0] < val):
		return np.nan
	else:
		while (b - a)/2.0 > tol:
			midpoint = (a + b)/2.0
			D,M = dense(R2,midpoint)
			if D < val: # Increasing but below 0 case
				b = midpoint
				mass = M
			else:
				a = midpoint
		return(midpoint,mass)



######################################################################
#%% CALCULATE BOUND AND UNBOUND PARTICLES
from numba import cuda

#a_device = cuda.to_device(a)
#b_device = cuda.to_device(b)
#out_device = cuda.device_array(shape=(row,col,), dtype=np.float64)

#@vectorize(['float64(float64, float64)'], target='cuda')
#def fast_pot_E(x,y):
	#return x**2 + y**2



#@guvectorize([(float64[:],float64, float64[:])], '(n,m),(n,k)->(n)',nopython=True)
#def net_vec(coord,res):
#	for i in prange(c):
#		res[i] = np.sum(coord - c[i])
#		res[i] = np.sum(c-c[i],axis=0)/r.size
		
	

@njit(parallel=True,fastmath=True)
def pot_E(c):
	potential = np.empty_like(c)
		#pointer = m**2*np.nanmean((c-c[i]),axis=0)
		#vec = pointer/np.linalg.norm(pointer)
		#Shifted points with CoM:  (CoM[:,np.newaxis,:]-c)	 
		#potential[i] = sum(coord_diff) / len(coord_diff)
		#potential[i] = np.sum((c-c[i]),axis=0)
		#print(i,end='\r')
	for i in prange(c.shape[0]):
		potential[i] = np.sum(c-c[i],axis=0)/c.shape[0]
	#potential = net_vec(c)
	potential = Grav_const*m**2*np.sqrt(potential[:,0]**2+potential[:,1]**2+potential[:,2]**2)
	return potential

pot_E.parallel_diagnostics(level=4)

	

def kin_E(m_part,v_part,v_com):
	return 0.5*m_part*np.linalg.norm(v_part-v_com,axis=1)**2


def bound_checker(m,v,c,r,v_com):#,masses,velo,Cosmology):
	#implement cosmology later! Returns True/False array for each particle
	EK = kin_E(m,v,v_com)
	EP = pot_E(c)
	print('EK',EK,'EP',EP)
	bound_states = np.where(EK<EP,True,False)
	#bound_states = np.ones_like(radii, dtype=bool)
	return np.array(bound_states)


######################################################################
#%% Calculating Overdensity Radii & Plotting
rm_200m_list = []
rm_500c_list = []
timestep_list = []
filenamez = np.array([])
nframes = 1

for index, radii in enumerate(radius_particles):
	framenr = index
	timestep_list.append(times[index])
	rm_500c = bisec(R_min,R_max,radii,500*rho_c,1E-13)
	rm_500c_list.append(rm_500c)
	rm_200m = bisec(R_min,R_max,radii,200*rho_m,1E-13)
	rm_200m_list.append(rm_200m)
	nan_mask = np.array(~np.isnan(radii))	#Ignoring the nans!
	bound = bound_checker(masses[index][:][nan_mask],new_velos_list[index][:][nan_mask],new_coords_list[index][:][nan_mask]-CoM[index],radii[nan_mask],CoV[index]) 
	print('\t Bound %= ',bound.sum()/bound.size*100,'%',end='\r')
	filenamez = np.append(filenamez,np.array([output_foldername+'images/'+filename+str(framenr)]*nframes))
	dens = density(rsp)
	c = densities[index] > 0
	dens = np.compress(c, densities[index])
	rs_i = np.compress(c, rs)

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
		x=new_coords_list[index,:,0][bound],
		y=new_coords_list[index,:,1][bound],
		z=new_coords_list[index,:,2][bound],
	#Unbound:
		x1=new_coords_list[index,:,0][bound != True],
		y1=new_coords_list[index,:,1][bound != True],
		z1=new_coords_list[index,:,2][bound != True],
	#Sphereical Overdensity:
		centre=CoM[index,:],
		#(shift*2,shift*2),
		radi_200m=rm_200m[0],
		m200m=rm_200m[1],
		radi_500c=rm_500c[0],
		m500c=rm_500c[1],
	#Plotting Info:
		lims=[-0.1*np.array(f[0]["Header"].attrs["BoxSize"])[0], 1.1*np.array(f[0]["Header"].attrs["BoxSize"])[0]],
		folder=output_foldername+'images/',
		name=filename+str(framenr),
		title='Spatial distribution of halo particles',
		rs=rs_i,
		dens=dens,
		Rlims=(Rmin,Rmax),
		Denslims=(Dmin,Dmax),
		dens_model=[])

rm_500c_array = np.array(rm_500c_list)
rm_200m_array = np.array(rm_200m_list)
timestep = np.array(timestep_list).ravel()

plot_evo(timestep,rm_500c_array[:,0],r'R$_{500,c}$',rm_200m_array[:,0],r'R$_{200,m}$','Evolution of the overdensity radii','Time [Gyr]',r'Radius [kpc h$^{-1}$]',output_foldername)
plot_evo(timestep,rm_500c_array[:,1],r'M$_{500,c}$',rm_200m_array[:,1],r'M$_{200,m}$','Evolution of the overdensity masses','Time [Gyr]',r'Mass [10$^{10}$ M$_{\odot}$ h$^{-1}$]',output_foldername)


#sys.exit()

######################################################################
#% Build GIF        
print('creating gif density profile\n')
with imageio.get_writer(output_foldername+'Overdensity region plot.gif', mode='I') as writer:
	for fn in filenamez:
		print('Loading: '+str(fn),end='\r')
		ims = imageio.v2.imread(str(fn)+'.png')
		writer.append_data(ims)
		os.remove(fn+'.png')
print('gif complete\n')
print('Done!')


