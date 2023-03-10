from numpy import random, cos, sin, sqrt, log
import sys
import os
sys.path.insert(0,'ProfileFiles/')
import time,math
import pandas as pd
from ICs_NFW import *
from ICs_NFWX import *
from ICs_Hernquist import *
from ICs_King import *
from ICs_Einasto import *


'''''''''''''''''''''''''''''''''''''''''''''''''''''''''
'''''''''''''''''''''''''''''''''''''''''''''''''''''''''
'''''''''''''''''''''''''''''''''''''''''''''''''''''''''
				READ PARAMETER FILE
'''''''''''''''''''''''''''''''''''''''''''''''''''''''''
'''''''''''''''''''''''''''''''''''''''''''''''''''''''''
'''''''''''''''''''''''''''''''''''''''''''''''''''''''''
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
    m = d['m']
    G = d['G']
    params = d
    return model,n,m,G,params

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
    while not good_n:
    	if trun and model == 'NFW':
    		n=int(1.43*target_n**1.008+100)
    	elif model == 'Hernquist':
    		n = int(target_n)
    
    	#Find R
	rand = random.rand(n)
	initialguess = eval(model+'_RootGuess')(params)
	res = optimize.root(eval(model+'_MD'),initialguess*rand,args=(rand,params),tol=1e-4,method='krylov')
	R = res.x

	#Find V
	P = eval(model+'_P')(R,params) #potential
	E = pick_E(model,P,params) #energy
	V = sqrt(2*(P-E)) #velocity magnitude

	#Convert dimensionless R,V to r,v
	r,v = eval(model+'_dimens')(R,V,n,m,G,params)
	
	if r.size >= target_n:
		good_n = True
    
    print('Used Particles: '+str(100-(float((r.size-target_n))/target_n*100))+'%')    
    #x,y,z components of positions and velocities
    sort_radius = np.argsort(r)
    #r = np.sort(r)
    
    r = r[sort_radius]
    v = v[sort_radius]
    
    x,y,z = isotrop(r[:target_n])
    vx,vy,vz = isotrop(v[:target_n])
    
    return x, y, z, vx, vy, vz

	
def pick_E(model,P,params):
#########################################################
##Given points with potentials P, and a distribution 
##function, DF, picks random energies from the distribution
#########################################################


    #Generate distribution function
    E = np.linspace(0.0,max(P),1e3)
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
model,n,m,G,params = read_paramfile(paramfile)
eval(model+'_Initialize')(params)

#Run Program
x, y, z ,vx, vy, vz = ICs(model,n,m,G,params)

#Write File
filename = sys.argv[2]
if filename[-4:] == '.txt':
    filename = filename[:-4]+'.csv'

try:
    os.remove(filename)
except (OSError,TypeError,ValueError) as e:
    pass

d = {'x': x, 'y': y,'z': z, 'vx': vx,'vy': vy, 'vz': vz}
df = pd.DataFrame(data=d)

#df = df*1000
posi = df[['x','y','z']]
min_pos = posi[['x','y','z']].min().min()
max_pos = posi[['x','y','z']].max().max()
max_box = round((abs(min_pos)+abs(max_pos))*2,6)
centre_box = max_box/2

df[['x','y','z']] = df[['x','y','z']] + centre_box

#print(df.shape)
df.to_csv(filename,sep=',')

with open('Params.txt', 'w') as f:
    f.write(str(max_box)+" "+str(m))
f.close()



