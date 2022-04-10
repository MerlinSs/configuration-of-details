#!/usr/bin/python
# -*- coding: utf-8 -*-
import glob
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib.font_manager import FontProperties

import os

######################## Plotting ####################################
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
mpl.rcParams['xtick.direction'] = 'in'
mpl.rcParams['ytick.direction'] = 'in'
mpl.rcParams['xtick.top'] = True
mpl.rcParams['ytick.right'] = True
nrows = 2
ncols = 1
fig, axs = plt.subplots(nrows, ncols, sharex = True, figsize=(12,8))

txt_size = 22
    
###### Subroutine to read the correlators and the n_k0 #############################
def read_E_0(filename):

	E_0 = np.nan
    
	try:
		f_in = open(filename, 'r')
		#print (filename)
		str_to_find = 'Ground'

		for line in f_in:
			if  str_to_find in line:
				string_split = line.split()
				try:
					E_0 = float(string_split[-1]) 
				except ValueError:
					pass
				break
		f_in.close()       
	except IOError:
		print (" %s not found"  %(filename))
         
	return E_0

def getD(phi_x_values, V_0):

	M_E_0 = np.full((len(V_0), len(phi_x_values)), np.nan)
	L_E_0 = np.full(len(phi_x_values), np.nan)
	D = np.full(len(V_0), np.nan)
	
	for j in range(len(V_0)):

		for i_phi_x, phi_x in enumerate(phi_x_values):

			filename = path + "stripe_lx%dly%dP%d_%dhcbost%.1fV%.2fV_0_%.5fphi_x%.5fphi_y%.5f_nx%dny%d.out" %(lx, ly, P, nhcbos, t, V, V_0[j], phi_x, phi_y, nx, ny)
			L_E_0[i_phi_x] = read_E_0(filename)
		
		M_E_0[j, :] = L_E_0
		
		fits = np.polyfit(phi_x_values, M_E_0[j], 2)
		D[j] = 2*fits[0]
	
	#print(M_E_0)
	
	return D
	
def read_correlators(filename, nev):

    E_alpha = np.full(nev, np.nan)
    real_c_20 = np.full(nev, np.nan)
    real_c_21 = np.full(nev, np.nan)
    real_n_k0 = np.full(nev, np.nan)
    
    try:
        f_in = open(filename, 'r')
        
        
    
        str_to_find = ' Eigenvalues \t c_20 \t c_21 \t n_k0'
        
        i_nev = 0
        for line in f_in:
            if str_to_find in line:
                    
                str_to_stop = '##############################'
                for line_read in f_in:
                    #print (line_read)
                    if str_to_stop in line_read:
                        break
                        
                    string_split = line_read.split()
                        
                    if not string_split:  ## This is to check if arrived in the empty line
                        break
                        
                    if i_nev == nev:
                        break
                            
                    try:
                            
                        E_alpha[i_nev] = float(string_split[0])
                        real_c_20[i_nev] = float(string_split[1])
                        real_c_21[i_nev] = float(string_split[6])
                        real_n_k0[i_nev] = float(string_split[11])

                        imag_c_20 = float(string_split[3])
                        imag_c_21 = float(string_split[8])
                        imag_n_k0 = float(string_split[13])
                            
                        if (np.fabs(imag_c_20) > 1e-12):
                            raise("imag_c_20 is too large")
                        if (np.fabs(imag_c_21) > 1e-12):
                            raise("imag_c_21 is too large")
                        if (np.fabs(imag_n_k0) > 1e-12):
                            raise("imag_n_k0 is too large")
                        i_nev+=1
                            
                    except ValueError:
                        pass
                break
        
            
        f_in.close()
        
    except IOError:
        print ('File not found - %s' %(filename))

    return E_alpha, real_c_20, real_c_21, real_n_k0

def read_densities(filename, P):

    dens = np.full(P, np.nan)
    
    try:
        f_in = open(filename, 'r')
        
        str_to_find = 'reference site\t<n_i>'
        
        idx = 0
        for line in f_in:
            if str_to_find in line:
                
                for line_read in f_in:
                    #print (line_read)
                        
                    string_split = line_read.split()
                        
                    if not string_split:  ## This is to check if arrived in the empty line
                        break
                        
                    if idx == P:
                        break
                            
                    try:
                            
                        #print (float(string_split[1]))
                        dens[idx] = float(string_split[1])
                        idx+=1
                        
                    except ValueError:
                        pass
                break
        
            
        f_in.close()
        
    except IOError:
        print ('File not found - %s' %(filename))
        
    return dens
    
def read_c_01_correlator(filename):
    
    c_01 = np.nan
    
    reference_site = 3 
    dx = 0
    dy = 1
    
    try:
        f_in = open(filename, 'r')
        
        str_to_find = 'reference site\tdx\tdy\tninj'
        for line in f_in:
            if str_to_find in line:
                
                str_to_find_2 = '%lu\t%lu\t%lu\t' %(reference_site, dx, dy) 
                for line_read in f_in:
                    if str_to_find_2 in line_read:
                        string_split = line_read.split()
                        
                        try:
                            c_01 = float(string_split[-1])
                        except ValueError:
                            pass
                        break
                break
    except IOError:
        print ('File not found - %s' %(filename))
        
    return c_01

################ Parameters ###############################################################
#### Lattice geometry #####################################################################
lx = 8
ly = 4
P = 4

#### Filling ##############################################################################
nhcbos = 16

#### Hamiltonian parameters ###############################################################
t = 1.00000
phi_x = 0.00000
phi_x_values = 0.000
#phi_x_values = np.arange(0, np.pi - 0.001, 0.7853)
phi_y = 0.00000
V = 4.00000000
list_of_V_0s = np.unique(np.concatenate(([0.1, 0.01], np.arange(0, 1, 0.1), np.arange(0, 11, 1), np.arange(10,30,2), np.arange(10,101,10))))

### Momentum sectors ######################################################################
nx = 0
ny = 0

kx = 2.0*nx*np.pi/(lx)
ky = 2.0*ny*np.pi/(ly)

#### Number of eigenvalues computed #######################################################
nev = 4

path = '/home/jh/got/hc_%dbosons/' %(nhcbos)
#path = '/home/jh/codes/hc_%dbosons/' %(nhcbos)

### Formatting lists #####
list_of_colors = ['dodgerblue', 'crimson', 'forestgreen', 'purple']
list_of_markers = ['o', 's', '^', 'v']
list_of_labels = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h']

### Allocating room for the data ##################
E_alpha =  np.full([len(list_of_V_0s), nev], np.nan)
c_20 = np.full([len(list_of_V_0s), nev], np.nan)
c_21 = np.full([len(list_of_V_0s), nev], np.nan)
n_k0 = np.full([len(list_of_V_0s), nev], np.nan)
c_01 = np.full(len(list_of_V_0s), np.nan)
c_01_test = np.full(len(list_of_V_0s), np.nan)
dens = np.full([len(list_of_V_0s), P], np.nan)
#D = getD(phi_x_values, list_of_V_0s)
#print(D)

for i_V_0, V_0 in enumerate(list_of_V_0s):
    
    ## Read correlators ###
    filename = path + "stripe_lx%luly%luP%lu_%luhcbost%.1fV%.2fV_0_%.5fphi_x%.5fphi_y%.5f_nx%luny%lu.out" %(lx,ly,P,nhcbos,t,V,V_0,phi_x,phi_y,nx,ny)
    E_alpha[i_V_0, :], c_20[i_V_0, :], c_21[i_V_0, :], n_k0[i_V_0, :] = read_correlators(filename, nev)
    
    ## Read inhomogeneous densities ###
    filename = path + "stripe_hc_bosons_lx%luly%luP%lu_%luhcbost%.1fV%.2fV_0_%.5fphi_x%.5fphi_y%.5f_kx%.5fky%.5fE_0.out" %(lx,ly,P,nhcbos,t,V,V_0,phi_x,phi_y,kx,ky)
    dens[i_V_0, :] = read_densities(filename, P)
    c_01[i_V_0] = read_c_01_correlator(filename)

#####################################origenal print(include D)#########################
#axs[0].plot(list_of_V_0s, c_20[:, 0] - dens[:, 1]*dens[:, 1], color = list_of_colors[0], marker = list_of_markers[0], linestyle = '', markersize = 8, fillstyle = 'none', label = r'$c(2,0)$')
#axs[0].plot(list_of_V_0s, c_21[:, 0] - dens[:, 1]*dens[:, 1], color = list_of_colors[1], marker = list_of_markers[1], linestyle = '', markersize = 8, fillstyle = 'none', label = r'$c(2,1)$')
#twin_1, = axs[1].plot(list_of_V_0s, n_k0[:, 0], color = list_of_colors[2], marker = list_of_markers[2], linestyle = '', markersize = 8, fillstyle = 'none', label = r'$n_{k=0}$')

#twin2 = axs[1].twinx()
#twin_2, = twin2.plot(list_of_V_0s, D, color = 'purple', linestyle = '-.', marker = 'v', markersize = 8, fillstyle = 'none', label = r'$D$')

#axs[1].legend(handles = [twin_1, twin_2], prop={'size': 22})

#axs[0].axhline(0.0, ls = '--', lw = 1, color = 'dimgrey')     

##### My initial convention ############################
#axs[0].plot(list_of_V_0s, c_20[:, 0] - dens[:, 1]*dens[:, 1], color = list_of_colors[0], marker = list_of_markers[0], linestyle = '', markersize = 8, fillstyle = 'none', label = r'$c(2,0)$')
#axs[0].plot(list_of_V_0s, c_21[:, 0] - dens[:, 1]*dens[:, 1], color = list_of_colors[1], marker = list_of_markers[1], linestyle = '', markersize = 8, fillstyle = 'none', label = r'$c(2,1)$')

##### Another convention #################################
axs[0].plot(list_of_V_0s, c_20[:, 0] - 0.5*(dens[:, 1] + dens[:, 3]) + 0.25, color = list_of_colors[0], marker = list_of_markers[0], linestyle = '', markersize = 8, fillstyle = 'none', label = r'$c_3(2,0)$')
axs[0].plot(list_of_V_0s, c_21[:, 0] - 0.5*(dens[:, 1] + dens[:, 3]) + 0.25, color = list_of_colors[1], marker = list_of_markers[1], linestyle = '', markersize = 8, fillstyle = 'none', label = r'$c_3(2,1)$')
axs[0].plot(list_of_V_0s, c_01[:] - 0.5*(dens[:, 3] + dens[:, 3]) + 0.25, color = list_of_colors[2], marker = list_of_markers[2], linestyle = '', markersize = 8, fillstyle = 'none', label = r'$c_3(0,1)$')

axs[1].plot(list_of_V_0s, n_k0[:, 0], color = list_of_colors[3], marker = list_of_markers[3], linestyle = '', markersize = 8, fillstyle = 'none')

axs[0].axhline(0.0, ls = '--', lw = 1, color = 'dimgrey')

####### fine tuning plot ###################################    
axs[nrows-1].set_xlabel(r'$V_0/t$', fontsize = txt_size)

for i in range(nrows):
    axs[i].tick_params(axis='both', which='major', labelsize= 22, length=10)
    # Initialize minor ticks
    axs[i].minorticks_on()
        
axs[0].set_ylabel(r'$c({\bf r})$', fontsize = txt_size)
axs[1].set_ylabel(r'$n_{{\bf k}=0}$', fontsize = txt_size)
#twin2.set_ylabel(r'$D = \frac{\partial E^2_0}{\partial^2 \phi_x}$', fontsize = txt_size)

#### making labels for each panel ###########
axs[0].text(0.1, 0.1, r'$({\bf %s})$' %(list_of_labels[0]), transform = axs[0].transAxes, fontsize = txt_size)
axs[1].text(0.1, 0.1, r'$({\bf %s})$' %(list_of_labels[1]), transform = axs[1].transAxes, fontsize = txt_size)
#axs[0, 1].text(0.1, 0.1, r'$({\bf %s})$' %(list_of_labels[2]), transform = axs[0, 1].transAxes, fontsize = txt_size)
#axs[1, 1].text(0.1, 0.1, r'$({\bf %s})$' %(list_of_labels[3]), transform = axs[1, 1].transAxes, fontsize = txt_size)

fig.suptitle(r'$%d \times %d\ \ {\cal P} = %d\ \ \ $' %(lx, ly, P) + r'$V/t = %.1f \ \ n_{\rm hcbos} = %d$' %(V, nhcbos), fontsize = txt_size)

#### Making legend ##########################################
axs[0].legend(loc='best', numpoints = 1, fontsize = txt_size)
              #ncol=4, mode="expand", borderaxespad=0., fontsize = 0.8*txt_size)

#### Whole figure finie tuning ##############################
#plt.subplots_adjust(hspace = 0.0, wspace = 0.4)

############### Outputting figure ###########################
#output_path='/home/jh/图片/research_print/'
#output_path='/home/jh/codes/figure/'
#file_name = 'n=%lu,E_0,ref_s=%d,ref_s_test=%d.png' %(nhcbos, ref_s, ref_s_test) 
#file_name = os.path.join(output_path, file_name)
#plt.savefig(file_name, dpi='figure', orientation='portrait', format="png", bbox_inches=None)

plt.show()
