# -*- coding: utf-8 -*-
# IMPORT THE IMPORTANT PACKAGES FOR PYTHON
from sklearn.linear_model import LinearRegression
from scipy.optimize import curve_fit
import numpy as np
import pandas as pd
import h5py
import matplotlib.pyplot as plt
import scipy.optimize as sco
import sys
import os 
np.random.seed(42069)

#PLOT PARAMS

# Set the font size for axis labels and tick labels
plt.rcParams['font.size'] = 20

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
plt.rcParams['figure.dpi'] = 200

# Set the default file format to PDF
plt.rcParams['savefig.format'] = 'pdf'
plt.rcParams.update({'axes.titlesize': 'large'})

# Set the default dpi of figures to 150
#plt.rcParams['figure.dpi'] = 200



# LOAD ALL THE FILES FOR PLOTTING:
model = str(sys.argv[1])
files = np.array(os.listdir('Results'))

am_files = files[np.flatnonzero(np.core.defchararray.find(files,'Angular_momentum')!=-1)]
hv_files = files[np.flatnonzero(np.core.defchararray.find(files,'Halo_vision')!=-1)]

filef = open(f"../Final_Results/{model}/Out_final.txt", "w")

# Define power law function
def power_law(x, a, b):
    return a * x ** b

# define the uncertainties
def calc_sig(y_lower_error,y_upper_error):
	sigma = np.zeros_like(y_lower_error)
	sigma[y_upper_error > y_lower_error] = y_upper_error[y_upper_error > y_lower_error]
	sigma[y_upper_error <= y_lower_error] = y_lower_error[y_upper_error <= y_lower_error]
	return sigma

#Plotting function(s)
def init_plot(params):
	title,xlab,ylab = params
	fig, ax = plt.subplots()
	ax.set_title(title,fontsize=18)
	ax.set_xlabel(xlab)
	ax.set_ylabel(ylab)
	return ax
	
def plotter(axs,data,col,alp,lab):
	xdata,ydata = zip(*list(data))
	axs.plot(xdata,ydata,'o-',c=col,alpha=alp,label=lab)

def fin_plot(axs,folder,title,lims=False):
	axs.legend(loc='lower left',fancybox=True, shadow=True,fontsize=16)
	axs.tick_params(axis='both', which='minor', length=3, labelcolor='black',bottom=True, top=True, left=True, right=True)
	axs.tick_params(axis='both', which='major', length=5, labelcolor='black',bottom=True, top=True, left=True, right=True)
	if lims:
		axs.set_xlim(left=-0.2)
		axs.set_ylim(bottom=0.45,top=1.04)
	fig = axs.get_figure()
	fig.set_tight_layout(True)
	fig.savefig(folder+title+'.png',bbox_inches='tight')
	plt.close(fig)
	#fig.show()
	#plt.pause(0.00001)

colors = ['blue', 'orange', 'green']

# PLOTTING PHASE:
xaxlabs 	= ["Time [Gyr]",'Time [Gyr]','Time [Gyr]',"Time [Gyr]","Time [Gyr]"]
plot_titles = ["Ratio Evolution of the overdensity masses","Evolution of the overdensity masses","Evolution of the overdensity radii","Energies","Angular Momentum"]
yaxlabs 	= ["Ratio",r'Mass $\times$ 10$^{10}$ $h^{-1}$ [M$_{\odot}$]',r'Radius $h^{-1}$ [kpc]',"log Energy [..]",r"Angular momentum [$10^{10}$ M$_{\odot}$ kpc$^{2}$ s$^{-1}$]"]

for massa in [12,13,14]:
	am_ma = am_files[np.flatnonzero(np.core.defchararray.find(am_files,f'{massa}')!=-1)]
	hv_ma = hv_files[np.flatnonzero(np.core.defchararray.find(hv_files,f'{massa}')!=-1)]
	mass_ejected	= []
	ejection_time	= []
	extra_mass_lost_200m	= []
	extra_mass_lost_500c	= []
	ejection_err_200m 	= []
	ejection_err_500c 	= []
	mass500_ejected		= []
	for ejet in ['-4-','-5.5-','-7-']:
		am_ej = am_ma[np.flatnonzero(np.core.defchararray.find(am_ma,f'{ejet}')!=-1)]
		hv_ej = hv_ma[np.flatnonzero(np.core.defchararray.find(hv_ma,f'{ejet}')!=-1)]
		#print(am_ej)
		axes = []			
		solm = '\odot'
		ejet = ejet.replace('-','')
		params = list(zip([fr'M$_{{200m}}$ = $10^{{{massa}}}$ $h^{{-1}}$ M$_{solm}$ | Event: {ejet} Gyr' for pt in plot_titles],xaxlabs,yaxlabs))
		alp = 1
		labels_00 = set()
		labels_01 = set()
		lab200m = fr'200$_m$'
		lab500c = fr'500$_c$'
		labx = 'X'
		laby = 'Y'
		labz = 'Z'
		labt = 'Total'
		inits = True
		for ejem in range(0,11):#['-2-','-5-','-10-']:
			ejem = '-'+str(ejem)+'-'
			am_ps = am_ej[np.flatnonzero(np.core.defchararray.find(am_ej,f'{ejem}{ejet}')!=-1)]
			hv_ps = hv_ej[np.flatnonzero(np.core.defchararray.find(hv_ej,f'{ejem}{ejet}')!=-1)]
			#for parsel in ['-F']:
			parsel = '-F'
			#i = str(am_ps[np.flatnonzero(np.core.defchararray.find(am_ps,f'{parsel}')!=-1)].item())
			#j = str(hv_ps[np.flatnonzero(np.core.defchararray.find(hv_ps,f'{parsel}')!=-1)].item())
			parsel = parsel.replace('-','')
			ejem = ejem.replace('-','')
			am_df = pd.read_csv(f'Results/{am_ps.item()}')
			hv_df = pd.read_csv(f'Results/{hv_ps.item()}')
			data_200m_r = zip(hv_df.Times,np.array(hv_df.m200m)/np.array(hv_df.m200m)[0])
			data_500c_r = zip(hv_df.Times,np.array(hv_df.m500c)/np.array(hv_df.m500c)[0])
			if len(axes) < 4:
				axes.append(init_plot(params[0]))
			plotter(axes[0],data_200m_r,'r',alp,lab=lab200m if lab200m not in labels_00 else '')
			plotter(axes[0],data_500c_r,'b',alp,lab=lab500c if lab500c not in labels_00 else '')
			data_200m = zip(hv_df.Times,hv_df.m200m)
			data_500c = zip(hv_df.Times,hv_df.m500c)
			if len(axes) < 4:
				axes.append(init_plot(params[1]))
			plotter(axes[1],data_200m,'r',alp,lab=lab200m if lab200m not in labels_00 else '')
			plotter(axes[1],data_500c,'b',alp,lab=lab500c if lab500c not in labels_00 else '')
			data_200m = zip(hv_df.Times,hv_df.radi_200m)
			data_500c = zip(hv_df.Times,hv_df.radi_500c)
			if len(axes) < 4:
				axes.append(init_plot(params[2]))
			plotter(axes[2],data_200m,'r',alp,lab=lab200m if lab200m not in labels_00 else '')
			plotter(axes[2],data_500c,'b',alp,lab=lab500c if lab500c not in labels_00 else '')
			#Ten = zip(hv_df.Times,hv_df.Total_Energy)
			#Pen = zip(hv_df.Times,hv_df.Potential_Energy)
			#Ken = zip(hv_df.Times,hv_df.Kinetic_Energy)
			#if len(axes) < 4:
			#	axes.append(init_plot(params[3]))
			#plotter(axes[2],Ten,'purple',alp,lab=fr'Total E |{parsel} {ejem}%')
			#plotter(axes[2],Pen,'b',alp,lab=fr'PE |{parsel} {ejem}%')
			#plotter(axes[2],Ken,'r',alp,lab=fr'KE |{parsel} {ejem}%')
			amx = zip(am_df.Time,np.log10(am_df.Total_x_angmom))
			amy = zip(am_df.Time,np.log10(am_df.Total_y_angmom))
			amz = zip(am_df.Time,np.log10(am_df.Total_z_angmom))
			amt = zip(am_df.Time,np.log10(am_df.Total_angmom))
			if len(axes) < 4:
				axes.append(init_plot(params[3]))
			plotter(axes[3],amx,'b',alp,lab=labx if labx not in labels_01 else '')
			plotter(axes[3],amy,'y',alp,lab=laby if labx not in labels_01 else '')
			plotter(axes[3],amz,'g',alp,lab=labz if labx not in labels_01 else '')
			plotter(axes[3],amt,'r',alp,lab=labt if labx not in labels_01 else '')	
			alp -= 0.085
			data_200m_r = np.array(hv_df.m200m)/np.array(hv_df.m200m)[0]
			data_500c_r = np.array(hv_df.m500c)/np.array(hv_df.m500c)[0]
			ejection_time.append(float(ejet))
			mass_ejected.append(int(ejem))
			EJT_200m = round(data_200m_r.shape[0]*float(ejet) /14)
			EJT_relax = 30
			if inits:
				null_meet = data_200m_r
				null_meetr = data_500c_r
				inits = False
			before_200m = np.median(data_200m_r[:EJT_200m]-null_meet[:EJT_200m]+1)
			after_200m =  np.median(data_200m_r[EJT_200m+EJT_relax:]-null_meet[EJT_200m+EJT_relax:]+1)
			mloss_200m =  data_200m_r[EJT_200m-1]-data_200m_r[EJT_200m]
			#print(mloss_200m)
			#print(before_200m,after_200m,mloss_200m,ejem,before_200m-after_200m-mloss_200m,np.percentile(before_200m-data_200m_r[EJT_200m+EJT_relax:]-mloss_200m,(16,50,84)))
			extra_mass_lost_200m.append(before_200m-after_200m-mloss_200m)
			ejection_err_200m.append(np.percentile(before_200m-data_200m_r[EJT_200m+EJT_relax:]-null_meet[EJT_200m+EJT_relax:]+1-mloss_200m,(16,50,84)))
			EJT_500c = round(data_500c_r.shape[0]*float(ejet) /14)
			before_500c = np.median(data_500c_r[:EJT_500c]-null_meetr[:EJT_500c]+1)
			after_500c =  np.median(data_500c_r[EJT_500c+EJT_relax:]-null_meetr[EJT_500c+EJT_relax:]+1)
			mloss_500c =  data_500c_r[EJT_500c-1]-data_500c_r[EJT_500c]
			extra_mass_lost_500c.append(before_500c-after_500c-mloss_500c)
			ejection_err_500c.append(np.percentile(before_500c-data_500c_r[EJT_500c+EJT_relax:]-null_meetr[EJT_500c+EJT_relax:]+1-mloss_500c,(16,50,84)))
			mass500_ejected.append(np.round(np.abs((hv_df.m200m[EJT_200m-1]*int(ejem)/100)/(hv_df.m500c[EJT_500c-1])*100),2))
			if lab200m not in labels_00:
				labels_00.add(lab200m)
				labels_00.add(lab500c)
			if labx not in labels_01:
				labels_01.add(labx)
				labels_01.add(laby)
				labels_01.add(labz)
				labels_01.add(labt)
		for i,axs in enumerate(axes):
			fin_plot(axs,
				folder = f"../Final_Results/{model}/",
				title = fr'{model} {plot_titles[i]} ({massa}_{parsel}_{ejet})',
				)#lims=True)
	#plt.hist(before-mloss-data_200m_r[EJT_200m:],bins=50)
	#plt.hist(before-mloss-data_200m_r[EJT_200m:],bins=50)
	#plt.show()
	print('plotting')
	ejection_time = np.array(ejection_time) 
	mass_ejected = np.array(mass_ejected) 
	extra_mass_lost = np.array(extra_mass_lost_200m)	
	ejection_err = np.array(ejection_err_200m) #do i need to subtract the error??
	#print(ejection_err)
	figu, ax = plt.subplots()
	# Create a set to keep track of the unique labels
	labels_added = []
	tmass = str(massa)
	#ax.set_title(r'DM Halo Extra R$_{200m}$ OD-Mass Loss from Ejection (M$_{200m}$ = '+r'10$^{' + str(massa) + r'}$'+r' M$_{\odot}$)')
	ax.set_xlabel(r'Halo Mass Ejected [%]')
	ax.set_ylabel(r'Extra Overdensity Massloss [%]')
	x  = np.unique(mass_ejected)
	q = x.size
	xx = x[1:]
	y1 = (extra_mass_lost[1:q])*100
	y2 = (extra_mass_lost[q+1:2*q])*100
	y3 = (extra_mass_lost[2*q+1:])*100
	y_err1 = (ejection_err[1:q])*100
	y_err2 = (ejection_err[q+1:2*q])*100
	y_err3 = (ejection_err[2*q+1:])*100

	y_err1s = calc_sig(y_err1[:,1]-y_err1[:,0],y_err1[:,2]-y_err1[:,1])
	y_err2s = calc_sig(y_err2[:,1]-y_err2[:,0],y_err2[:,2]-y_err2[:,1])
	y_err3s = calc_sig(y_err3[:,1]-y_err3[:,0],y_err3[:,2]-y_err3[:,1])
	
	y_err1 = np.array(list(zip(y_err1[:,1]-y_err1[:,0],y_err1[:,2]-y_err1[:,1]))).T
	y_err2 = np.array(list(zip(y_err2[:,1]-y_err2[:,0],y_err2[:,2]-y_err2[:,1]))).T
	y_err3 = np.array(list(zip(y_err3[:,1]-y_err3[:,0],y_err3[:,2]-y_err3[:,1]))).T


# print the fit parameters and their uncertainties


	# Plot the data with error bars and create the legend
	pwlb = []
	for i, (y, y_err, y_errs, label) in enumerate(zip([y1, y2, y3], [y_err1, y_err2, y_err3], [y_err1s, y_err2s, y_err3s], [f'{ejection_time[0]}', f'{ejection_time[q]}', f'{ejection_time[2*q]}'])):
		popt, pcov = curve_fit(power_law, xx, y)#, sigma=y_errs)
		x = np.append(0,xx)
		y = np.append(0,y)
		y_err = np.concatenate((np.array([0,0]).reshape(-1,1),y_err),axis=1)
		ax.errorbar(x-0.1+i*0.1, y, yerr=y_err, color=colors[i], fmt='o', capsize=4, label=f'{label} Gyr' if label not in labels_added else '')
		# Generate points for fitted power law curve
		x_fit = np.linspace(1E-10, 10.5+0.1-i*0.1, 100)
		y_fit = power_law(x_fit, *popt)
		# Calculate error regions		
		popt_u, pcov_u = curve_fit(power_law, xx, y[1:], sigma=y_err[1,1:].flatten())#,absolute_sigma=True)
		y_fit_upper = power_law(x_fit, *(popt + np.sqrt(np.diag(pcov))))
		popt_l, pcov_l = curve_fit(power_law, xx, y[1:], sigma=y_err[0,1:].flatten())#,absolute_sigma=True)
		y_fit_lower = power_law(x_fit, *(popt - np.sqrt(np.diag(pcov))))
		print("a = {:.3f} +/- {:.3f}".format(popt[0], np.sqrt(pcov[0, 0])))
		print("b = {:.3f} +/- {:.3f}".format(popt[1], np.sqrt(pcov[1, 1])))
		filef.write('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n')
		filef.write(f'{massa} results (200m), {label}\n')
		filef.write("a = {:.3f} +/- {:.3f}\n".format(popt[0], np.sqrt(pcov[0, 0])))
		filef.write("b = {:.3f} +/- {:.3f}\n".format(popt[1], np.sqrt(pcov[1, 1])))
		bb = f'{popt[1]:.3f}'
		pwlb.append(fr'y = {popt[0]:.3f} x'+r'$^{' + str(bb) + r'}$')
		# Plot the fitted curve and error regions
		ax.plot(x_fit-0.1+i*0.1, y_fit,'--', color=colors[i], label=f'{label} Gyr')
		label = f'{label} Gyr'
		#plt.fill_between(x_fit-0.1+i*0.1, y_fit_lower, y_fit, color=colors[i], alpha=0.1,label='Fit error')
		#plt.fill_between(x_fit-0.1+i*0.1, y_fit, y_fit_upper, color=colors[i], alpha=0.1)
		if label not in labels_added:
			labels_added.append(label)
	ax.set_xlim(left=0.0)
	ax.set_ylim(bottom=0.0)
	ax.tick_params(axis='both', which='minor', length=3, labelcolor='black',bottom=True, top=True, left=True, right=True)
	ax.tick_params(axis='both', which='major', length=5, labelcolor='black',bottom=True, top=True, left=True, right=True)
	eb = ax.get_legend_handles_labels()[0][3:]
	#print(pwlb)
	labels = ["Ejection Time"] + labels_added + ["Power Law Fit"] + pwlb
	ph = [plt.plot([],marker="", ls="")[0]]*2
	handles = ph[:1] + eb + ph[1:] + ax.lines[3::4]
	leg = plt.legend(handles, labels, ncol=2, loc='upper left',fancybox=True, shadow=True,fontsize=16)
	ax.set_title(fr'M$_{{200m}}$ (M$_{{200m}}$ = $10^{{{massa}}}$ $h^{{-1}}$ M$_{solm}$)')
	figu.set_tight_layout(True)
	figu.savefig(f"../Final_Results/{model}/{model}: {massa} results (200m).png",bbox_inches='tight')
	plt.close(figu)

	extra_mass_lost = np.array(extra_mass_lost_500c)
	ejection_err = np.array(ejection_err_500c)
	figu, ax = plt.subplots()
	# Create a set to keep track of the unique labels
	labels_added = []
	tmass = str(massa)
	#ax.set_title(r'DM Halo Extra R$_{500c}$ OD-Mass Loss from Ejection (M$_{200m}$ = '+r'10$^{' + str(massa) + r'}$'+r' M$_{\odot}$)')
	ax.set_xlabel(r'Halo Mass Ejected [%]')
	ax.set_ylabel(r'Extra Overdensity Massloss [%]')
	x  = np.unique(mass_ejected)
	q = x.size
	x = np.array(mass500_ejected)

	xx1 = (x[1:q])
	xx2 = (x[q+1:2*q])
	xx3 = (x[2*q+1:])

	y1 = (extra_mass_lost[1:q])*100
	y2 = (extra_mass_lost[q+1:2*q])*100
	y3 = (extra_mass_lost[2*q+1:])*100

	y_err1 = (ejection_err[1:q])*100
	y_err2 = (ejection_err[q+1:2*q])*100
	y_err3 = (ejection_err[2*q+1:])*100

	y_err1s = calc_sig(y_err1[:,1]-y_err1[:,0],y_err1[:,2]-y_err1[:,1])
	y_err2s = calc_sig(y_err2[:,1]-y_err2[:,0],y_err2[:,2]-y_err2[:,1])
	y_err3s = calc_sig(y_err3[:,1]-y_err3[:,0],y_err3[:,2]-y_err3[:,1])
	
	y_err1 = np.array(list(zip(y_err1[:,1]-y_err1[:,0],y_err1[:,2]-y_err1[:,1]))).T
	y_err2 = np.array(list(zip(y_err2[:,1]-y_err2[:,0],y_err2[:,2]-y_err2[:,1]))).T
	y_err3 = np.array(list(zip(y_err3[:,1]-y_err3[:,0],y_err3[:,2]-y_err3[:,1]))).T

	
	# Plot the data with error bars and create the legend
	pwlb = []
	for i, (xx, y, y_err, y_errs, label) in enumerate(zip([xx1, xx2, xx3], [y1, y2, y3], [y_err1, y_err2, y_err3], [y_err1s, y_err2s, y_err3s], [f'{ejection_time[0]}', f'{ejection_time[q]}', f'{ejection_time[2*q]}'])):
		popt, pcov = curve_fit(power_law, xx, y)#, sigma=y_errs)
		x = np.append(0,xx)
		y = np.append(0,y)
		y_err = np.concatenate((np.array([0,0]).reshape(-1,1),y_err),axis=1)
		ax.errorbar(x-0.2+i*0.2, y, yerr=y_err, color=colors[i], fmt='o', capsize=4, label=f'{label} Gyr' if label not in labels_added else '')
		# Generate points for fitted power law curve
		x_fit = np.linspace(1E-10, np.ceil(np.max(xx+0.2)), 100)
		y_fit = power_law(x_fit, *popt)
		# Calculate error regions		
		popt_u, pcov_u = curve_fit(power_law, xx, y[1:], sigma=y_err[1,1:].flatten())#,absolute_sigma=True)
		y_fit_upper = power_law(x_fit, *(popt + np.sqrt(np.diag(pcov))))
		popt_l, pcov_l = curve_fit(power_law, xx, y[1:], sigma=y_err[0,1:].flatten())#,absolute_sigma=True)
		y_fit_lower = power_law(x_fit, *(popt - np.sqrt(np.diag(pcov))))
		print("a = {:.3f} +/- {:.3f}".format(popt[0], np.sqrt(pcov[0, 0])))
		print("b = {:.3f} +/- {:.3f}".format(popt[1], np.sqrt(pcov[1, 1])))
		filef.write('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n')
		filef.write(f'{massa} results (500c), {label}\n')
		filef.write("a = {:.3f} +/- {:.3f}\n".format(popt[0], np.sqrt(pcov[0, 0])))
		filef.write("b = {:.3f} +/- {:.3f}\n".format(popt[1], np.sqrt(pcov[1, 1])))
		bb = f'{popt[1]:.3f}'
		pwlb.append(fr'y = {popt[0]:.3f} x'+r'$^{' + str(bb) + r'}$')
		# Plot the fitted curve and error regions
		plt.plot(x_fit-0.2+i*0.2, y_fit,'--', color=colors[i], label=f'{label} Gyr')
		#plt.fill_between(x_fit-0.1+i*0.1, y_fit_lower, y_fit, color=colors[i], alpha=0.1,label='Fit error')
		#plt.fill_between(x_fit-0.1+i*0.1, y_fit, y_fit_upper, color=colors[i], alpha=0.1)
		if label not in labels_added:
			labels_added.append(label)
	ax.set_xlim(left=0.0)
	ax.set_ylim(bottom=0.0)
	ax.tick_params(axis='both', which='minor', length=3, labelcolor='black',bottom=True, top=True, left=True, right=True)
	ax.tick_params(axis='both', which='major', length=5, labelcolor='black',bottom=True, top=True, left=True, right=True)
	h, l = ax.get_legend_handles_labels()
	eb = ax.get_legend_handles_labels()[0][3:]
	handles = ax.lines[1::4] + ax.lines[3::4]
	labels = ["Ejection Time"] + labels_added + ["Power Law Fit"] + pwlb
	ph = [plt.plot([],marker="", ls="")[0]]*2
	handles = ph[:1] + eb + ph[1:] + ax.lines[3::4]
	leg = plt.legend(handles, labels, ncol=2, loc='upper left',fancybox=True, shadow=True,fontsize=16)
	ax.set_title(fr'M$_{{500c}}$ (M$_{{200m}}$ = $10^{{{massa}}}$ $h^{{-1}}$ M$_{solm}$)')
	figu.set_tight_layout(True)
	figu.savefig(f"../Final_Results/{model}/{model}: {massa} results (500c).png",bbox_inches='tight')
	plt.close(figu)

filef.close()

