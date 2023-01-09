# -*- coding: utf-8 -*-
################################################################################
#
# Run with: python profiles_gif.py Figures/ NFW
#
################################################################################

from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import numpy as np
import h5py as h5
import imageio
import sys
import os

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
params = d

data_foldername = 'SwiftRUN_'+model+'/'
output_foldername = str(sys.argv[1])


#%%
if not os.path.exists('Figures/images'):
	os.mkdir('Figures/images')


Nmax = 100
steps = 1
n_frames = 1

#%% Parameters
combine = True #input('Combine the Density & Mass Profiles in 1 subplot? (Y/N)') 

iterarray = np.arange(0, Nmax + 1, steps)
f = [h5.File(data_foldername+"output_%04d.hdf5" % i, "r") for i in iterarray]
coordinates_lists = [i["PartType1"]["Coordinates"] for i in f]
masses = [i["PartType1"]["Masses"] for i in f]
times = [i["Header"].attrs["Time"][0] * 
          i["Units"].attrs["Unit time in cgs (U_t)"][0] / 31557600.0e9 for i in f]

shift = f[0]["Header"].attrs["BoxSize"] / 2
nbsamples = 250
Rmin = 1
Rmax = np.linalg.norm(np.array(f[0]["Header"].attrs["BoxSize"]),axis=-1)
G = 4.299581e04 #??
figsize = 10

xmin = 10**.8
xmax = f[0]["Header"].attrs["BoxSize"][0]*10**.5
ymin0 = 10**-13.2
ymax0 = 10**-3.8
ymin1 = np.array(masses[0]).min()*10**-.2
ymax1 = np.array(masses[0]).sum()*10**.2


min_d,max_d = 0,0
min_m,max_m = 0,0
rsp = np.logspace(np.log10(Rmin), np.log10(Rmax), nbsamples)
#rsp = np.linspace((Rmin), (Rmax), nbsamples)

M = np.sum(masses[0])
filenames1, filenames2 = [],[]

time_array = np.zeros(len(iterarray))

#radius_particles = np.linalg.norm(np.array(coordinates_lists)-shift,axis=-1)
min_radius = 0.001 #radius_particles.min() 
max_radius = Rmax #np.linalg.norm(shift) #radius_particles.max()

def overdens(m,r):
	rho_c = 9*10**-18 #kg/k*m**3 #3H**2/8*np.pi*Grav_const #2.64419905 *10**23 #kg/pc**3
	mass_bins = np.array([np.array(masses)[radius_particles<r].sum() for r in rsp])
	densities = mass_bins/(4/3*np.pi*rsp**3)
	densities[densities==0] = np.nan
	r_overdens = rsp[np.argwhere(densities<=9*10**-18)][0]
	# np.array[masses[q]/(4/3*radius_particles[q]**3*np.pi) for q in range(Nmax)])

def r_vir(H,G,factor,crit=True):
	if crit:
		rho=...
	else:
		rho=...
	M_vir=4/3*np.pi*r
	
	return virial_radius,virial_mass
	

def mass_ins(R):
	return ((r < R) * mass).sum()
	
def density(R):
	return np.diff(mass_ins_vect(R)) / np.diff(R) / (4.0 * np.pi * R[1:] ** 2)

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


#dens model	= eval("model_"+model)
#mass_cum_model	= dens model



for i in iterarray:
	f = h5.File(data_foldername+"output_%04d.hdf5" % i, "r")
#	f = h5.File("output_%04d.hdf5" % i, "r")   

	time_array[int(i / steps)] = f["Header"].attrs["Time"]
	boxsize = 	f["Header"].attrs["BoxSize"] / 2
	particles = 	f["PartType1"]
	coordinates = 	particles["Coordinates"][:, :]
	velocities = 	particles["Velocities"][:, :]
	masses = 	particles["Masses"][:]

	pos = np.array(f["DMParticles"]["Coordinates"]) - shift

	time = (
	f["Header"].attrs["Time"][0]
	* f["Units"].attrs["Unit time in cgs (U_t)"][0]
	/ 31557600.0e9
	)

	mass = np.array(f["DMParticles"]["Masses"])
	r = np.linalg.norm(np.array(pos),axis=-1)

	# Methods to compute density profile

	mass_ins_vect = np.vectorize(mass_ins)
	#m_array_cum = np.array([mass[r<k].sum() for k in arr])

	dens = density(rsp)
	rs = rsp[1:]

	# remove empty bins
	c = dens > 0
	dens = np.compress(c, dens)
	rs = np.compress(c, rs)

	mass_ = np.diff(mass_ins_vect(rsp))
	c1 = mass_ > 0
	mass_ = np.compress(c1, mass_)

	mass_cum = mass_ins_vect(rsp)

	filename1 = f'images/frame_{i}_D.png'
	filename2 = f'images/frame_{i}_M.png'
	#print(filename1)
	for frame in range(n_frames):
		filenames1.append(filename1)
		filenames2.append(filename2)	
	
	
# Plotting the figures
	
	if combine:
		fig, ax = plt.subplots(2, sharex=True,figsize=(1.2 * figsize, figsize))#, constrained_layout=True)
		fig.suptitle(r"$t=$ {:.3f} Gyr".format(time), fontsize=16)
		#ax[0].set(xlim=(xmin,xmax),ylim=(ymin0,ymax0))
		#ax[1].set(xlim=(xmin,xmax),ylim=(ymin1,ymax1))
		ax[0].set_title(
		    r"Density Profile ($M = {:.1e}$ M$_{{\odot}}$)".format(
		        M * 1e10
		    )
		)
		ax[1].set_title(
		    r"Mass Profile ($M = {:.1e}$ M$_{{\odot}}$)".format(
		        M * 1e10
		    )
		)
		ax[0].set_ylabel(r"$\rho(r)$ [$M_{\odot}$ kpc$^{-3}$]")
		ax[1].set_ylabel(r"M [$M_{\odot}$]")
		ax[0].plot(rs, dens, label=r"Simulation")
		#ax[0].plot(rs, dens_model, label=r"Model",c='r',ls='-')		
		ax[1].plot(rs, mass_, label=r"Radial Mass")
		ax[1].plot(rsp, mass_cum, label='Cumulative mass')
		#ax[1].plot(rsp, mass_cum_model, label='Cumulative mass (model)',c='r',ls='-')

#		parameters, covariance = curve_fit(plummer_analytical, rs, dens)
#		a = parameters[0]

#		
#		ax[0].plot(rsp, plummer_analytical(rsp,a),ls='--',linewidth=3,alpha=.6, c="black", label="Analytical Model")
#		ax[1].plot(rsp, m_model(rsp,a),ls='--',linewidth=3,alpha=.6, c="black", label="Analytical Cumulative Model")
#		ax[1].plot(rsp[1:], np.diff(m_model(rsp,a)),ls='--',linewidth=3,alpha=.6, c="red", label="Analytical Shell Model")
		ax[1].set_xlabel("r [kpc]")
		ax[0].legend()
		ax[1].legend()
		ax[0].loglog()
		ax[1].loglog()
		plt.tight_layout()
		#print(i,output_foldername+filename1)
		plt.savefig(output_foldername+filename1)
		plt.close()
	    
	else:
		#Dens plot
		fig, ax = plt.subplots(1,figsize=(1.2 * figsize, figsize))
		ax.plot(rs, dens, label=r"$t=$ {:.3f} Gyr".format(time))
	    
		parameters, covariance = curve_fit(plummer_analytical, rs, dens)
		a = parameters[0]
		
		ax.plot(rsp, plummer_analytical(rsp,a),ls='--',linewidth=3,alpha=.6, c="black", label="Analytical Model")
		ax.set_xlabel("r [kpc]")
		ax.legend()
		ax.loglog()
		ax.set_title(
		    r"Plummer Density Profile: $a = {:.1e}$ kpc, $M = {:.1e}$ M$_{{\odot}}$".format(
		        a, M * 1e10
		    )
		)
		plt.tight_layout()
		ax.set_ylabel(r"$\rho(r)$ [$M_{\odot}$ kpc$^{-3}$]")
		plt.savefig(f'images/frame_{i}_D.png')
		plt.close()
		
		#Mass plot
		fig1, ax1 = plt.subplots(1,figsize=(1.2 * figsize, figsize))
		ax1.plot(rs, dens, label=r"$t=$ {:.3f} Gyr".format(time))
	    
		parameters, covariance = curve_fit(plummer_analytical, rs, dens)
		a = parameters[0]
		
		ax.plot(rsp, plummer_analytical(rsp,a),ls='--',linewidth=3,alpha=.6, c="black", label="Analytical Model")
		ax.set_xlabel("r [kpc]")
		ax.legend()
		ax.loglog()
		ax.set_title(
		    r"Plummer Mass Profile: $a = {:.1e}$ kpc, $M = {:.1e}$ M$_{{\odot}}$".format(
		        a, M * 1e10
		    )
		)
		plt.tight_layout()
		ax.set_ylabel(r"$\rho(r)$ [$M_{\odot}$ kpc$^{-3}$]")
		plt.savefig(output_foldername+filename2)
		plt.close()


#% Build GIF        
print('creating gif density profile\n')
#print(filenames1)
with imageio.get_writer(output_foldername+'dens_profile (Model: '+model+').gif', mode='I') as writer:
    for filename in filenames1:
        image = imageio.v2.imread(output_foldername+filename)
        writer.append_data(image)

print('gif complete\n')
print('Removing Images\n')

# Remove files
#for filename in set(filenames1):
#    os.remove(output_foldername+filename)
print('done #1')

# Build GIF for mass
if not combine:        
    print('creating gif mass profile\n')
    with imageio.get_writer(output_foldername+'mass_profile (Model: '+model+').gif', mode='I') as writer:
        for filename in filenames2:
            image = imageio.v2.imread(output_foldername+filename)
            writer.append_data(image)
    
    print('gif complete\n')
    print('Removing Images\n')
    
    # Remove files
#    for filename in set(filenames2):
#        os.remove(output_foldername+filename)
    print('done #2')
