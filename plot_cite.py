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
    
###################################Function############ #############################
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

def read_correlator(filename, ref_site):
    
    M = np.full((4, 8), np.nan)
    
    reference_site = ref_site
    dx = 0
    dy = 0
    
    try:
        f_in = open(filename, 'r')
        
        str_to_find = 'reference site\tdx\tdy\tninj'
        
        for line in f_in:
            if str_to_find in line:
                for dy in np.arange(4):
                    for dx in np.arange(8):
                        str_to_find_2 = '%lu\t%lu\t%lu\t' %(reference_site, dx, dy) 
                
                        for line_read in f_in:
                            if str_to_find_2 in line_read:
                                string_split = line_read.split()
                        
                                try:
                                    print(float(string_split[-1]))
                                    M[dx, dy] = float(string_split[-1])
                                except ValueError:
                                    pass
    except IOError:
        print ('File not found - %s' %(filename))
        
    return M


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
phi_y = 0.00000
V = 1.00000000
V_0 = 0.0000 
ref_site = 2
### Momentum sectors ######################################################################
nx = 0
ny = 0

kx = 2.0*nx*np.pi/(lx)
ky = 2.0*ny*np.pi/(ly)

#### Number of eigenvalues computed #######################################################
nev = 4

path = './hc_16bosons/'

### Formatting lists #####
list_of_colors = ['dodgerblue', 'crimson', 'forestgreen', 'purple']
list_of_markers = ['o', 's', '^', 'v']
list_of_labels = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h']

### Allocating room for the data ##################
dens = np.full(P, np.nan)
M = np.full((4, 8), np.nan)
filename = path + "stripe_hc_bosons_lx%luly%luP%lu_%luhcbost%.1fV%.2fV_0_%.5fphi_x%.5fphi_y%.5f_kx%.5fky%.5fE_0.out" %(lx,ly,P,nhcbos,t,V,V_0,phi_x,phi_y,kx,ky)
M = read_correlator(filename, ref_site)
dens = read_densities(filename, P)
print(M, dens)
