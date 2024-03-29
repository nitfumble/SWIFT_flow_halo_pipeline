import numpy as np
from numpy import sqrt, pi, log, inf
from scipy import integrate, optimize, special

'''''''''''''''''''''''''''''''''''''''''''''''''''''''''
'''''''''''''''''''''''''''''''''''''''''''''''''''''''''
'''''''''''''''''''''''''''''''''''''''''''''''''''''''''
				NFWX FUNCTIONS
'''''''''''''''''''''''''''''''''''''''''''''''''''''''''
'''''''''''''''''''''''''''''''''''''''''''''''''''''''''
'''''''''''''''''''''''''''''''''''''''''''''''''''''''''

def NFWX_Initialize(params):
#########################################################
##Initialize
#########################################################
    #Check decay parameter OK
    d = params['d']
    c = params['r_vir']/params['r_sX']
    criteria = (1+c)*(1+c)*((1+c)*log(1+c)-c)/(1.0+3.0*c)/(2.0*(1+c)*log(1+c)-c)
    if d<criteria:
        print('Warning: Your decay parameter is too small. May return unphysical ICs.\n')
    return


def NFWX_RootGuess(params):
#########################################################
##Returns initial guess for root finder method
#########################################################

    r_s = params['r_sX']
    r_vir = params['r_vir']
    c = r_vir/r_s
    
    return c



def NFWX_Mass(R,c,d):
#########################################################
## Mass within radius R for NFW profile with exponential cutoff
########################################################

    #Make sure R is an array
    if not isinstance(R, np.ndarray):
        R = np.array([R])

    #Mass interior to viral radius
    m = np.minimum(R,c)
    m = np.maximum(m,0) #so doesn't throw error in root-finder
    M = log(m+1)-m/(m+1)

    #Add in mass past virial radius (if R>c)
    epsilon = -(1+3*c)/(1.0+c) + c/d
    subgamma = special.gammaincc(epsilon+3,c/d) - special.gammaincc(epsilon+3,R[R>c]/d)
    subgamma = special.gamma(epsilon+3)*subgamma #since pythons gamma function normalized
    M[R>c] += (d/(1+c))**2 * (d/c)**(epsilon+1) * np.exp(c/d)*subgamma

    return M


def NFWX_MD(R,x,params):
#########################################################
##Cummulative mass distribution for NFWX profile
##Given random number x, root gives radius R
#########################################################

    r_s = params['r_sX']
    r_vir = params['r_vir']
    c = r_vir/r_s
    d = params['d']
    epsilon = -(1.0+3*c)/(1.0+c) + c/d
    Mtot = (d/(1.0+c))**2 * (d/c)**(epsilon+1) * np.exp(c/d)*special.gamma(epsilon+3)*special.gammaincc(epsilon+3,c/d)
    Mtot += log(c+1.0)-c/(c+1.0)

    MD = NFWX_Mass(R,c,d) - x*Mtot

    return MD

def NFWX_P(R,params):
#########################################################
##Potential in an NFWX profile
#########################################################
    #Make sure R is an array
    # print('R:',R)
    if not isinstance(R, np.ndarray):
        R = np.array([R])
    
    r_s = params['r_sX']
    r_vir = params['r_vir']
    c = r_vir/r_s
    d = params['d']

    #Infinite part of integral (maxR to inf)
    def integrand(x,c,d):
        return NFWX_Mass(x,c,d)/(x*x)
        
    Pend = integrate.quad(integrand,max(R),inf,args=(c,d))[0]

    #Integrate until maxR value
    x = np.linspace(min(R),max(R),int(1e4))
    # print('x:',x)
    y = NFWX_Mass(x,c,d)/(x*x)
    Ps = integrate.cumtrapz(y,x,initial=0) #gives integral from rmin to r
    Ps = Ps[-1]-Ps #int_rmin^rmax - int_rmin^r = int_r^rmax
    # print('Ps:',Ps)
    P = np.interp(R.ravel(),x.ravel(),Ps.ravel())+Pend
    
    return P

def NFWX_DF(E,params):
#########################################################
##Distribution function at relative energy E
##For NFWX profile
#########################################################
    Emin = min(E/10.0,1e-1) #pick how far to integrate out to 

    r_s = params['r_sX']
    r_vir = params['r_vir']
    c = r_vir/r_s
    d = params['d']
    
    #Find rmin
    P0 = NFWX_P(1E-5,params)
    def myfunc(R,E):
        if R<0: #so root finder doesn't find negative R
            P = P0 - R
        else:
            P = NFWX_P(R,params)
        return E -P
    res = optimize.root(myfunc,(1-E)/E,args=(E,),tol=1e-4)
    Rmin = res.x
    Rmax = optimize.root(myfunc,1.0/Emin,args=(Emin,),tol=1e-4) #radii at minimum E
    Rmax = Rmax.x
    
    def integrand(R,params):
        epsilon = -(1+3*c)/(1+c) + c/d
        L = NFWX_Mass(R,c,d)
        
        p = 1.0/(R*(1+R)**2)
        p[R>c] = 1.0/c/(1.0+c)**2 *(R[R>c]/c)**epsilon * np.exp(-(R[R>c]-c)/d)
        
        f = (R/(L*(1+R)))**2 *p* (6*L + p*R*(1+R)*(3*R+1))
        f[R>c] = p[R>c]/L[R>c]/d/d *(d*d*epsilon*(epsilon -1) - 2*d*epsilon*R[R>c] + R[R>c]*R[R>c] +d*(epsilon*d - R[R>c])*(2 - p[R>c]*R[R>c]**3/L[R>c]))
        return f

    R = np.linspace(Rmin,Rmax,int(1e4)).reshape(-1)
    # print('linR:',R)
    P = NFWX_P(R,params)
    
    R = R[P<E]; P = P[P<E]

    f = integrand(R,params)/sqrt(E-P)
    F = integrate.simps(f,R)
    

    return F

    
def NFWX_dimens(R,V,n,m,G,params):
#########################################################
##Turns radius and velocity into dimension-full quantities
#########################################################
    r_s = params['r_sX']
    r_vir = params['r_vir']
    c = r_vir/r_s
    d = params['d']
    epsilon = -(1+3*c)/(1+c) + c/d
    
    Ltot = (d/(1+c))**2 * (d/c)**(epsilon+1) * np.exp(c/d)*special.gamma(epsilon+3)*special.gammaincc(epsilon+3,c/d)
    Ltot += log(c+1)-c/(c+1.0)
    
    p0 = n*m/(4.0*pi*r_s**3*Ltot)
                 
    r = R*r_s
    v = V*sqrt(4*pi*G*p0*r_s**2)

    return r,v

