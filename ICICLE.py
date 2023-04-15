from numpy import random, cos, sin, sqrt, log
import sys
sys.path.insert(0,'../ProfileFiles/')
#sys.path.insert(0,'ProfileFiles/')
import pandas as pd
import os
from ICs_NFW import *
from ICs_NFWX import *
from ICs_Hernquist import *
from ICs_King import *
from ICs_Einasto import *
import matplotlib.pyplot as plt

h = 0.6766
rhom = 8.634164473613977e-09
rhoc = 2.7753662724570803e-08

#PLOT PARAMS

# Set the font size for axis labels and tick labels
plt.rcParams['font.size'] = 18

# Set the font family for all text in the plot
plt.rcParams['font.family'] = 'serif'

# Set the figure size to 6 x 4 inches
plt.rcParams['figure.figsize'] = [8, 6]

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
plt.rcParams['figure.dpi'] = 300

# Set the default file format to PDF
plt.rcParams['savefig.format'] = 'pdf'

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
'''''''''''''''''''''''''''''''''''''''''''''''''''''''''
'''''''''''''''''''''''''''''''''''''''''''''''''''''''''
'''''''''''''''''''''''''''''''''''''''''''''''''''''''''
				READ PARAMETER FILE
'''''''''''''''''''''''''''''''''''''''''''''''''''''''''
'''''''''''''''''''''''''''''''''''''''''''''''''''''''''
'''''''''''''''''''''''''''''''''''''''''''''''''''''''''


def radii(rads,denz):
    a= rads.min()
    b= rads.max()
    global m
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
    print(midpoint, m * partii.size)
    return midpoint

def new_n(n,m):
    if m < 0.01:
        new_n_part = 145000 #360 26
    elif m < 0.1:
        new_n_part = 151700 #775 70
    else:
        new_n_part = 167600 #1670 220
    return new_n_part

def read_paramfile(paramfile):
#########################################################
##Read in Param File
#########################################################
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
    n = int(d['n'])
    trun = d['truncate']
    m_OG = d['m']
    trun = d['truncate']
    if model == 'Hernquist':
        if m_OG <= 0.1:
            d['m'] = d['m']/1.105
        else:
            d['m'] = d['m']/1.13
    if trun and model == 'NFW':
        d['m'] = d['m']*(n/new_n(n,m_OG))
    if model == 'Einasto':
        d['m'] = d['m']/1.5 
    if model == 'NFWX':
        d['m'] = d['m']/1.01    
    
    rvir = d['r_vir']
    m = d['m']
    G = d['G']
    params = d
    return model,n,m,G,params,m_OG,rvir

'''''''''''''''''''''''''''''''''''''''''''''''''''''''''
'''''''''''''''''''''''''''''''''''''''''''''''''''''''''
'''''''''''''''''''''''''''''''''''''''''''''''''''''''''
				MAIN PART OF PROGRAM
'''''''''''''''''''''''''''''''''''''''''''''''''''''''''
'''''''''''''''''''''''''''''''''''''''''''''''''''''''''
'''''''''''''''''''''''''''''''''''''''''''''''''''''''''

def ICs(model,n,m,G,params):
#########################################################
##Outputs positions and velocities of n particles in a
##profile specifed by model
#########################################################
    good_n = False
    target_n = n
    trun = params['truncate']
    global m_OG
    while not good_n:
        if trun and model == 'NFW':
            n=new_n(target_n,m_OG)
        if model == 'Einasto':
            n=1.5*target_n  
        if model == 'Hernquist':
            if m_OG <= 0.01:
                n=1.105*target_n
            else:
                n=1.13*target_n
        if model == 'NFWX':
            n=1.01*n

        #Find R
        rand = random.rand(int(n)).astype('f8')
        initialguess = eval(model+'_RootGuess')(params)
        res = optimize.root(eval(model+'_MD'),initialguess*rand,args=(rand,params),tol=1e-4,method='krylov')
        R = res.x
        
        #Find V
        P = eval(model+'_P')(R,params) #potential
        E = pick_E(model,P,params) #energy
        V = sqrt(2*(P-E)) #velocity magnitude
        #Convert dimensionless R,V to r,v
        r,v = eval(model+'_dimens')(R,V,n,m,G,params)
        
        print('try! '+str(r.size))
        if r.size >= target_n:
            good_n = True

    print('Used Particles: '+str(100-(float((r.size-target_n))/target_n*100))+'%')    
    #x,y,z components of positions and velocities
    sort_radius = np.argsort(r)
    sort_velos = np.argsort(v)
    #sort_radius = np.random.shuffle(sort_radius)

    r = r[sort_radius]
    v = v[sort_radius]
    # r = r[sort_velos]
    # v = v[sort_velos]

    # x,y,z = isotrop(r[:])
    # vx,vy,vz = isotrop(v[:])
    x,y,z = isotrop(r[:target_n])
    vx,vy,vz = isotrop(v[:target_n])

    return x, y, z, vx, vy, vz

	
def pick_E(model,P,params):
#########################################################
##Given points with potentials P, and a distribution 
##function, DF, picks random energies from the distribution
#########################################################


    #Generate distribution function
    E = np.linspace(0.0,max(P),int(1e3))
    f = np.zeros(E.shape[0])
    for i in range(1,E.shape[0]):
        f[i] = eval(model+'_DF')(E[i],params)

    #Make sure no nans
    E = E[~np.isnan(f)]
    f = f[~np.isnan(f)]

    #Choose E from distribution function
    n = P.shape[0]
    rand = random.rand(n)
    ans = np.ones(n)
    for i in range(n):
        newE = E[E<P[i]] #newE goes to maxE(R)=P(R)
        F = integrate.cumtrapz(f[E<P[i]]*sqrt(P[i]-newE),newE, initial=0) #cummulative distribution
        newE = newE[~np.isnan(F)]; F = F[~np.isnan(F)];
        ans[i] = np.interp(rand[i]*F[-1], F, newE) #find at which energy F=rand
   
    return ans
	

def isotrop(m):
#########################################################
##Given n points with radius or velocity magnitude, m, 
##outputs components in 3D for an isotropic distribution
#########################################################
    
	#Generate random numbers
	n = m.shape[0]
	u = 2*random.rand(n)-1 # uniform in -1, 1
	theta = 2*pi*random.rand(n) # uniform in 0, 2*pi

	#Find x,y,z components
	x = m*sqrt(1-u*u)*cos(theta)
	y = m*sqrt(1-u*u)*sin(theta)
	z = m*u

	return x, y, z

'''''''''''''''''''''''''''''''''''''''''''''''''''''''''
'''''''''''''''''''''''''''''''''''''''''''''''''''''''''
'''''''''''''''''''''''''''''''''''''''''''''''''''''''''
				MAKE OUTPUT FILE
'''''''''''''''''''''''''''''''''''''''''''''''''''''''''
'''''''''''''''''''''''''''''''''''''''''''''''''''''''''
'''''''''''''''''''''''''''''''''''''''''''''''''''''''''


#Initialize
paramfile = sys.argv[1]
#paramfile = 'params.txt'

model,n,m,G,params,m_OG,rvir = read_paramfile(paramfile)
eval(model+'_Initialize')(params)

# import the builtin time module
import time

# Grab Currrent Time Before Running the Code
start = time.time()
print('loading...')
#Run Program
x, y, z ,vx, vy, vz = ICs(model,n,m,G,params)

#Write File
filename = sys.argv[2]
#filename =  'out.txt'
if filename[-4:] == '.txt':
    filename = filename[:-4]+'.csv'

try:
    os.remove(filename)
except (OSError,TypeError,ValueError) as e:
    pass
#print(x.dtype)
d = {'x': x, 'y': y,'z': z, 'vx': vx,'vy': vy, 'vz': vz}
df = pd.DataFrame(data=d)
#print(df.dtypes)
#df = df*1000
posi = df[['x','y','z']]
min_pos = posi[['x','y','z']].min().min()
max_pos = posi[['x','y','z']].max().max()
max_box = round((abs(min_pos)+abs(max_pos))*2,6)
centre_box = max_box/2

df[['x','y','z']] = df[['x','y','z']] + centre_box

#print(df.shape)
df.to_csv(filename,sep=',',index=False)

df['r'] = np.sqrt(x**2+y**2+z**2)
rads = df.r
print(rads.min(),rads.max())

m = m_OG
if model == 'NFW':
	m = m_OG * n/np.sum(df.r<rvir)/1.014
if model == 'Hernquist':
	m = m_OG * n/np.sum(df.r<rvir)/1.001
if model == 'King':
	m = m_OG * n/np.sum(df.r<radii(rads,rhom*200))
if model == 'Einasto':
	if m_OG == 0.001:
		m = m_OG * n/np.sum(df.r<rvir)/1.16
	if m_OG == 0.01:
		m = m_OG * n/np.sum(df.r<rvir)/1.2
	if m_OG == 0.1:
		m = m_OG * n/np.sum(df.r<rvir)/1.61


print(f'rvir = {radii(rads,rhom*200)}')
print(f'r< = {np.sum(rads<rvir)}')
print(f'm = {m}')
with open('Params.txt', 'w') as f:
    f.write(str(max_box)+" "+str(m))
f.close()
# Grab Currrent Time After Running the Code
end = time.time()

#Subtract Start Time from The End Time
total_time = end - start
print("\n"+ str(int(total_time/60))+'m '+str(round(total_time%60,2))+'s')


#%%

fig = plt.figure()
plt.plot(x,y,'.',alpha=0.2,label='xy') 
plt.plot(x,z,'.',alpha=0.2,label='xz')
plt.plot(y,z,'.',alpha=0.2,label='yz')
plt.xlabel(r'Dimension [kpc]')
plt.ylabel(r'Dimension [kpc]')
plt.legend(loc='best')
plt.tight_layout()
# plt.show()
plt.savefig(f'Figures/IC_Dimensions ({m_OG*n}).png')
plt.close(fig)

fig1 = plt.figure()
plt.xlabel(r'Radius [kpc]')
plt.ylabel(r'Counts [#]')
plt.hist(rads,bins=100)
# plt.legend(loc='best')
plt.tight_layout()
# plt.show()
plt.savefig(f'Figures/IC_Histogram ({m_OG*n}).png')
plt.close(fig1)


#%%

import numpy as np

def density_profile(r, m, dr):
    """
    Calculates the density profile of a set of shells with mass m and radii r.

    Parameters:
    r (array): Radii of the shells.
    m (float): Total mass within each shell.
    dr (float): Width of each shell.

    Returns:
    density (array): Density profile of the shells.
    """
    volume = 4.0/3 * np.pi * (r**3)  # Calculate volume of a shell with radius r and width dr
    # volume = 4.0 * np.pi * (r**2) * dr  # Calculate volume of a shell with radius r and width dr
    density = m / volume                # Calculate density as mass / volume

    return density


# raddd = np.arange(rads.min(),rads.max(),1)  # Radii of the shells
raddd = np.logspace(np.log10(rads.min()),np.log10(rads.max()),1000)  # Radii of the shells
aa = np.ones_like(raddd)
for xi,kk in enumerate(raddd):
    aa[xi] *= np.sum(rads<kk)*m

mass = np.diff(aa)  # Total mass within each shell
dr = 0.5                            # Width of each shell

# density = density_profile(raddd[1:], mass, dr) # Calculate density profile
density = density_profile(raddd, aa, dr) # Calculate density profile
maskr = density > 0
density = np.compress(maskr, density)
raddd = np.compress(maskr, raddd)
mean_density = np.ones_like(raddd)*(rhom*200)
mean_cdensity = np.ones_like(raddd)*(rhoc*500)
#print(raddd[np.argwhere(np.diff(np.sign(density - mean_density)))][0])
#print(raddd[np.argwhere(np.diff(np.sign(density - mean_cdensity)))][0])


maskr[::25] = maskr[::25] * 0
density = np.compress(~maskr, density)
raddd = np.compress(~maskr, raddd)
mean_density = np.compress(~maskr, mean_density)
mean_cdensity = np.compress(~maskr, mean_cdensity)


#%%
fig2 = plt.figure()
plt.xlabel(r'Radius [kpc]')
plt.ylabel(r'$\rho$(<r) [M$_{\odot}$ kpc$^{-3}$]')
plt.loglog(raddd,mean_cdensity,'--',c='g',label=r'$\rho_{500c}$')
plt.loglog(raddd,mean_density,'--',c='r',label=r'$\rho_{200m}$')
plt.loglog(raddd,density,'o-',c='black',label=r'$\rho$ $_{IC}$')
plt.text(x=1,y=mean_density[0]*1.3,s=r'R$_{200m}$ = '+f'{radii(rads,rhom*200):.2f} kpc',c='r',size=14)
plt.text(x=1,y=mean_cdensity[0]*1.3,s=r'R$_{500c}$ = '+f'{radii(rads,rhoc*500):.2f} kpc',c='g',size=14)
plt.legend(loc='upper right')
plt.tight_layout()
# plt.show()
plt.savefig(f'Figures/IC_Density Profile ({m_OG*n}).png')
plt.close(fig2)

# Writing to file
with open(f"Figures/Data ({m_OG*n}).txt", "w") as file1:
    # Writing data to a file
    file1.write(f'R200m = {radii(rads,rhom*200)}\n')
    file1.write(f'R500c = {radii(rads,rhoc*500)}\n')


