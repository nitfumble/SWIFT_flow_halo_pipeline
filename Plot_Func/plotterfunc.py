#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 20 15:54:38 2022
@author: bart
"""

import numpy as np
from matplotlib import pyplot as plt
from matplotlib import patches as pch
from matplotlib.ticker import MultipleLocator

#plt.rcParams.update({'axes.labelsize' : 20})
#plt.rcParams.update({'axes.linewidth' : 2})
#plt.rcParams.update({'xtick.labelsize' : 14})
#plt.rcParams.update({'ytick.labelsize' : 14})
#plt.rcParams.update({'xtick.major.size' : 5})
#plt.rcParams.update({'ytick.major.size' : 5})
#plt.rcParams.update({'xtick.minor.size' : 3})
#plt.rcParams.update({'ytick.minor.size' : 3})
#plt.rcParams.update({'xtick.major.width' : 1.6})
#plt.rcParams.update({'ytick.major.width' : 1.6})
#plt.rcParams.update({'xtick.minor.width' : 1.2})
#plt.rcParams.update({'ytick.minor.width' : 1.2})
#plt.rcParams.update({'legend.fontsize' : 14})

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def plot3dhex(x,y,z,Gsize = 50):
    fig4, axs = plt.subplots(2, 2,figsize=(8,6),dpi =400 ,sharex='col', sharey='row')
    fig4.delaxes(axs[0][1])
    # fig4.subplots_adjust(top=0.9)
    fig4.suptitle('Data Distribution', fontsize=16, y=1.02, x=0.6)

    hb4 = axs[0, 0].hexbin(x,y,
                     gridsize=Gsize, bins='log', cmap='plasma',mincnt=1)

    axs[1, 0].hexbin(x,z,
                     gridsize=Gsize, bins='log', cmap='plasma',mincnt=1)

    axs[1, 1].hexbin(y,z,
                     gridsize=Gsize, bins='log', cmap='plasma',mincnt=1)


    axs[0, 0].set(ylabel='y [Mpc]')
    axs[1, 1].set(xlabel='y [Mpc]')
    axs[1, 0].set(xlabel='x [Mpc]',ylabel='z [pc]')

    axs[1][1].xaxis.set_tick_params(which='both', labelbottom=True, labeltop=False)

    cbaxes = fig4.add_axes([1.03, 0.1, 0.03, 0.85]) 
    cb = fig4.colorbar(hb4, ax=axs,  cax = cbaxes)
    cb.set_label('N')
    fig4.tight_layout()
    plt.show()

def plot3d_scat(x,y,z,col='black',Gsize = 50):
    fig4, axs = plt.subplots(2, 2,figsize=(8,6),dpi =320 ,sharex='col', sharey='row')
    fig4.delaxes(axs[0][1])
    fig4.suptitle('Data Distribution', fontsize=16, y=1.02, x=0.6)

    axs[0, 0].plot(x,y,'.',
                     c=col,linewidth=1)

    axs[1, 0].plot(x,z,'.',
                     c=col,linewidth=1)

    axs[1, 1].plot(y,z,'.',
                     c=col,linewidth=1)


    axs[0, 0].set(ylabel='y [Mpc]')
    axs[1, 1].set(xlabel='y [Mpc]')
    axs[1, 0].set(xlabel='x [Mpc]',ylabel='z [pc]')

    #axs[0, 0].xaxis.set_minor_locator(MultipleLocator(2.5))
    #axs[1, 0].xaxis.set_minor_locator(MultipleLocator(2.5))
    #axs[1, 1].xaxis.set_minor_locator(MultipleLocator(25))
    #axs[1, 0].yaxis.set_minor_locator(MultipleLocator(25))

    #axs[0, 0].legend(loc='best')



    axs[1][1].xaxis.set_tick_params(which='both', labelbottom=True, labeltop=False)

    fig4.tight_layout()
    plt.show()

def plot3d_2scat(x,y,z,x1,y1,z1,centre,radi_200m,m200m,radi_500c,m500c,folder,name,lims,col='black',col1='orange',Gsize = 50):
    fig4, axs = plt.subplots(2, 2,figsize=(10,10) ,sharex='col', sharey='row')
    fig4.delaxes(axs[0][1])
    #fig4.subplots_adjust(top=0.9)
    fig4.suptitle('Halo Structure '+str(name), fontsize=16)#, y=1.02, x=0.6)
    circle0m = pch.Circle((centre[0],centre[1]), radi_200m, color='r',lw=1.2,fill=False,zorder=5)
    circle1m = pch.Circle((centre[0],centre[2]), radi_200m, color='r',lw=1.2, fill=False,zorder=5)
    circle2m = pch.Circle((centre[1],centre[2]), radi_200m, color='r',lw=1.2, fill=False,zorder=5)

    circle0c = pch.Circle((centre[0],centre[1]), radi_500c, color='b',lw=1.2,fill=False,zorder=5)
    circle1c = pch.Circle((centre[0],centre[2]), radi_500c, color='b',lw=1.2, fill=False,zorder=5)
    circle2c = pch.Circle((centre[1],centre[2]), radi_500c, color='b',lw=1.2, fill=False,zorder=5)
    
    axs[0, 0].plot(x,y,'.',
                     c=col,linewidth=0.1,alpha=0.5,label='Bound')
    axs[0, 0].plot(x1,y1,'.',
                     c=col1,linewidth=0.1,alpha=0.5,label='Unbound')
    axs[0, 0].plot(centre[0],centre[1],'*',
                     c='pink',ms=10.0,alpha=1,label='Centre',zorder=10)
                     
    #Objects for in the legend:
    circle1 = axs[0, 0].plot([],[],'o', markerfacecolor='none', ms=14,
                     c='r',linewidth=1,label=r'R$_{200,m}$')
    circle2 = axs[0, 0].plot([],[],'o', markerfacecolor='none', ms=14,
                     c='b',linewidth=1,label=r'R$_{500,c}$')
    
    
    axs[0, 0].add_patch(circle0m)
    axs[0, 0].add_patch(circle0c)
    
    axs[0, 0].text(x=0.01, y=0.97,s=r'M$_{200,m}$= '+str(np.round(m200m,1))+r'×10$^{10}$ M$_{\odot}$/h',fontsize=11,c='r', weight="bold", ha='left', va='center', transform=axs[0, 0].transAxes)
    axs[0, 0].text(x=0.99, y=0.97,s=r'M$_{500,c}$= '+str(np.round(m500c,1))+r'×10$^{10}$ M$_{\odot}$/h',fontsize=11,c='b', weight="bold", ha='right', va='center', transform=axs[0, 0].transAxes)


    axs[1, 0].plot(x,z,'.',
                     c=col,linewidth=0.1,alpha=0.5)
    axs[1, 0].plot(x1,z1,'.',
                     c=col1,linewidth=0.1,alpha=0.5)
    axs[1, 0].plot(centre[0],centre[2],'*',
                     c='pink',ms=10.0,alpha=1,label='Centre',zorder=10)
    axs[1, 0].add_patch(circle1m)
    axs[1, 0].add_patch(circle1c)


    axs[1, 1].plot(y,z,'.',
                     c=col,linewidth=0.1,alpha=0.5)
    axs[1, 1].plot(y1,z1,'.',
                     c=col1,linewidth=0.1,alpha=0.5)
    axs[1, 1].plot(centre[1],centre[2],'*',
                     c='pink',ms=10.0,alpha=1,label='Centre',zorder=10)
    axs[1, 1].add_patch(circle2m)
    axs[1, 1].add_patch(circle2c)

#    axs[0, 0].set_ylim(np.nanmin(y)-0.25,np.nanmax(y)+0.25)
#    axs[0, 0].set_xlim(np.nanmin(x)-0.25,np.nanmax(x)+0.25)
#    axs[1, 0].set_xlim(np.nanmin(x)-0.25,np.nanmax(x)+0.25)
#    axs[1, 0].set_ylim(np.nanmin(z)-0.25,np.nanmax(z)+0.25)
#    axs[1, 1].set_xlim(np.nanmin(y)-0.25,np.nanmax(y)+0.25)

    axs[0, 0].set_ylim(lims)
    axs[0, 0].set_xlim(lims)
    axs[1, 0].set_xlim(lims)
    axs[1, 0].set_ylim(lims)
    axs[1, 1].set_xlim(lims)

    axs[0, 0].set(ylabel='y [Mpc/h]')
    axs[1, 1].set(xlabel='y [Mpc/h]')
    axs[1, 0].set(xlabel='x [Mpc/h]',ylabel='z [Mpc/h]')

    axs[1][1].xaxis.set_tick_params(which='both', labelbottom=True, labeltop=False)
    
    axs[0, 0].set(adjustable='box', aspect=1)
    axs[1, 0].set(adjustable='box', aspect=1)
    axs[1, 1].set(adjustable='box', aspect=1)
    
    #handlebox.add_artist(circle0c)
    #handlebox.add_artist(circle0m)
    
    legend = axs[0, 0].legend(#bbox_to_anchor=(0.5, 1.05),
    loc='lower right',
    fancybox=True, shadow=True, borderpad=0.2,
    framealpha=0.4,facecolor='grey')
    
    legend.legendHandles[0].set_ms(10)
    legend.legendHandles[1].set_ms(10)
    legend.legendHandles[2].set_ms(10)

    legend.legendHandles[0].set_alpha(1)
    legend.legendHandles[1].set_alpha(1)
    legend.legendHandles[2].set_alpha(1)
    
    fig4.tight_layout()
    plt.savefig(folder+str(name)+'.png',dpi=300)
    #plt.show()
    #plt.pause(1E-10)
    #plt.clf()
    #plt.cla()
    plt.close(fig4)
    
    
    
def plot_evo(x,data1,dat1lab,data2,dat2lab,title,xlab,ylab,folder):
	fig = plt.figure(figsize=(10,8))
	plt.title(title)
	plt.plot(x,data1,'o-',alpha=0.4,c='b',label=dat1lab)
	plt.plot(x,data2,'o-',alpha=0.4,c='r',label=dat2lab)
	plt.legend(loc='best',fancybox=True, shadow=True)
	plt.xlabel(xlab)
	plt.ylabel(ylab)
	plt.tight_layout()
	plt.savefig(folder+title+'.png')
	#plt.show()
	#plt.pause(0.1)
	#plt.clf()
	#plt.cla()
	plt.close(fig)
	return



def plot3d_2scat_prof(x,y,z,x1,y1,z1,centre,radi_200m,m200m,radi_500c,m500c,folder,name,lims,title,
					rs,dens,Rlims,Denslims,dens_model=[],col='black',col1='orange',Gsize = 50):
	fig4, axs = plt.subplots(2, 2,figsize=(10,10))
	axs[1, 0].sharex(axs[0, 0])
	axs[1, 1].sharex(axs[0, 0])
	axs[1, 0].sharey(axs[0, 0])
	axs[1, 1].sharey(axs[0, 0])
	#fig4.delaxes(axs[0][1])
	#fig4.subplots_adjust(top=0.9)
	fig4.suptitle(title, fontsize=16)#, y=1.02, x=0.6)
	circle0m = pch.Circle((centre[0],centre[1]), radi_200m, color='r',lw=1.2,fill=False,zorder=5)
	circle1m = pch.Circle((centre[0],centre[2]), radi_200m, color='r',lw=1.2, fill=False,zorder=5)
	circle2m = pch.Circle((centre[1],centre[2]), radi_200m, color='r',lw=1.2, fill=False,zorder=5)

	circle0c = pch.Circle((centre[0],centre[1]), radi_500c, color='b',lw=1.2,fill=False,zorder=5)
	circle1c = pch.Circle((centre[0],centre[2]), radi_500c, color='b',lw=1.2, fill=False,zorder=5)
	circle2c = pch.Circle((centre[1],centre[2]), radi_500c, color='b',lw=1.2, fill=False,zorder=5)

	axs[0, 0].plot(x,y,'.',
		             c=col,linewidth=0.1,alpha=0.5,label='Bound')
	axs[0, 0].plot(x1,y1,'.',
		             c=col1,linewidth=0.1,alpha=0.5,label='Unbound')
	axs[0, 0].plot(centre[0],centre[1],'*',
		             c='pink',ms=10.0,alpha=1,label='Centre',zorder=10)
		             
	#Objects for in the legend:
	circle1 = axs[0, 0].plot([],[],'o', markerfacecolor='none', ms=14,
		             c='r',linewidth=1,label=r'R$_{200,m}$')
	circle2 = axs[0, 0].plot([],[],'o', markerfacecolor='none', ms=14,
		             c='b',linewidth=1,label=r'R$_{500,c}$')


	axs[0, 0].add_patch(circle0m)
	axs[0, 0].add_patch(circle0c)

	axs[0, 0].text(x=0.01, y=0.97,s=r'M$_{200,m}$= '+str(np.round(m200m,1))+r'×10$^{10}$ M$_{\odot}$/h',fontsize=11,c='r', weight="bold", ha='left', va='center', transform=axs[0, 0].transAxes)
	axs[0, 0].text(x=0.99, y=0.97,s=r'M$_{500,c}$= '+str(np.round(m500c,1))+r'×10$^{10}$ M$_{\odot}$/h',fontsize=11,c='b', weight="bold", ha='right', va='center', transform=axs[0, 0].transAxes)

	axs[0, 1].set_title(r"Density Profile")# ($M = {:.1e}$ M$_{{\odot}}$)".format(M * 1e10))

	axs[0, 1].set_ylabel(r"$\rho(r)$ [$M_{\odot}$ kpc$^{-3}$]")

	axs[0, 1].plot(rs, dens,c='black', label=r"Simulation")
	axs[0, 1].plot(rs, dens,'.',ms=5,alpha=0.8,c='black', label=r"Simulation")

	#ax[0].plot(rs, dens_model, label=r"Model",c='r',ls='-o')		

	axs[0, 1].set_xlabel("r [kpc]")
	axs[0, 1].set(xlim=(Rlims),ylim=(Denslims))
	axs[0, 1].loglog()


	axs[1, 0].plot(x,z,'.',
		         c=col,linewidth=0.1,alpha=0.5)
	axs[1, 0].plot(x1,z1,'.',
		         c=col1,linewidth=0.1,alpha=0.5)
	axs[1, 0].plot(centre[0],centre[2],'*',
		         c='pink',ms=10.0,alpha=1,label='Centre',zorder=10)
	axs[1, 0].add_patch(circle1m)
	axs[1, 0].add_patch(circle1c)


	axs[1, 1].plot(y,z,'.',
		         c=col,linewidth=0.1,alpha=0.5)
	axs[1, 1].plot(y1,z1,'.',
		         c=col1,linewidth=0.1,alpha=0.5)
	axs[1, 1].plot(centre[1],centre[2],'*',
		         c='pink',ms=10.0,alpha=1,label='Centre',zorder=10)
	axs[1, 1].add_patch(circle2m)
	axs[1, 1].add_patch(circle2c)

#    axs[0, 0].set_ylim(np.nanmin(y)-0.25,np.nanmax(y)+0.25)
#    axs[0, 0].set_xlim(np.nanmin(x)-0.25,np.nanmax(x)+0.25)
#    axs[1, 0].set_xlim(np.nanmin(x)-0.25,np.nanmax(x)+0.25)
#    axs[1, 0].set_ylim(np.nanmin(z)-0.25,np.nanmax(z)+0.25)
#    axs[1, 1].set_xlim(np.nanmin(y)-0.25,np.nanmax(y)+0.25)

	axs[0, 0].set_ylim(lims)
	axs[0, 0].set_xlim(lims)
	axs[1, 0].set_xlim(lims)
	axs[1, 0].set_ylim(lims)
	axs[1, 1].set_xlim(lims)

	axs[0, 0].set(ylabel='y [Mpc/h]')
	axs[1, 1].set(xlabel='y [Mpc/h]')
	axs[1, 0].set(xlabel='x [Mpc/h]',ylabel='z [Mpc/h]')

	axs[1][1].xaxis.set_tick_params(which='both', labelbottom=True, labeltop=False)

	axs[0, 0].set(adjustable='box', aspect=1)
	axs[1, 0].set(adjustable='box', aspect=1)
	axs[1, 1].set(adjustable='box', aspect=1)
    
    #handlebox.add_artist(circle0c)
    #handlebox.add_artist(circle0m)
    
	legend = axs[0, 0].legend(#bbox_to_anchor=(0.5, 1.05),
	loc='lower right',
	fancybox=True, shadow=True, borderpad=0.2,
	framealpha=0.4,facecolor='grey')

	legend.legendHandles[0].set_ms(10)
	legend.legendHandles[1].set_ms(10)
	legend.legendHandles[2].set_ms(10)

	legend.legendHandles[0].set_alpha(1)
	legend.legendHandles[1].set_alpha(1)
	legend.legendHandles[2].set_alpha(1)

	fig4.tight_layout()
	plt.savefig(folder+str(name)+'.png',dpi=300)
	#plt.show()
	#plt.pause(1E-10)
	#plt.clf()
	#plt.cla()
	plt.close(fig4)
    
    
    





