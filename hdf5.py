# -*- coding: utf-8 -*-
"""
Created on Fri Apr 30 02:25:22 2021

@author: Jorge
"""

from qutip import *
import numpy as np
import matplotlib.pyplot as plt
import scipy.linalg as la

import h5py 
import os 

import filename_generator


path_data_QT = "../results/data/QT/Matrix_hdf5/" 
if not os.path.exists(path_data_QT): 
    os.makedirs(path_data_QT)
    
filename_QT = path_data_QT+filename_generator.filename_gen("QT",N,k0,d)
f = h5py.File(filename_QT, 'w')
g1 = f.create_group('1')
#M_num = M.states[0].full()
#len_M_i,len_M_j = n_points_x,n_points_y
#M_set = g1.create_dataset('M',(len(t_list),len_M_i,len_M_j),dtype='complex64')
dset_P1 = g1.create_dataset('P1',(len_Delta,len_Omega),dtype = 'float64')
dset_P2 = g1.create_dataset('P2',(len_Delta,len_Omega),dtype = 'float64')
dset_M = g1.create_dataset('M',(len_Delta,len_Omega),dtype = 'float64')
dset_Delta = g1.create_dataset('Delta',(len_Delta,),dtype = 'float64')
dset_Omega = g1.create_dataset('Omega',(len_Omega,),dtype = 'float64')
#dset_t = g1.create_dataset('t',(len(t_list)),dtype = 'float64')
#dset_t[:] = np.array(t_list)[:]

for t in range(0,len(t_list)):
    dset_M[t,:,:] = M.states[j].full[:,:]




