from qutip import *
import numpy as np
import h5py 
import os 

import plot
import filename_generator


#General Parameters

Gamma = 1
omega_points = 1024
omega_list = np.linspace(-200,200,omega_points)
epsilon = 10**(-2)
Delta = 0
N = 2
S = 2
Omega = 30*Gamma
level = 2
n_hat = [0,1,0]
k0=2*np.pi
k_vec = [0,k0,0]
k = np.linalg.norm(k_vec)
d = 0.05/k
sensor_method = 1
norm = 0

if sensor_method:
    path_plot = "../plots/data/maps/Photon_Spectrum_WS/PSC_N=%d_b_th_2"%(N)
    if not os.path.exists(path_plot): 
        os.makedirs(path_plot)
else:
    path_plot = "../plots/data/maps/Photon_Spectrum/PSC_N=%d_b_th_2"%(N)
    if not os.path.exists(path_plot): 
        os.makedirs(path_plot)
                             
if sensor_method:    
    path_data = "../results/data/QT/Photon_Spectrum/Matrix_PSCT_WS_hdf5_N=%d_S=%d_Omega=%.2f_Delta=%.2f_epsilon=%.2f/"%(N,S,Omega,Delta,epsilon)   
else:
    path_data = "../results/data/QT/Photon_Spectrum/Matrix_PSCT_hdf5_N=%d_S=%d_Omega=%.2f_Delta=%.2f_epsilon=%.2f/"%(N,S,Omega,Delta,epsilon)   


filename = path_data+filename_generator.filename_gen_sensor("PSCT",N,S,k0,d,Omega,Delta)
    
f = h5py.File(filename,'r')    
y_list = f['/1/y_list_data']            

#Parameter for the populations
F = N + S

x = omega_list

if norm:
    y_norm_list = []
    for i in range(0,len(y_list)):
        y_norm_i = (np.abs(y_list[i]) - np.abs(np.min(y_list)))/(np.abs(np.max(y_list))-np.abs(np.min(y_list)))
        y_norm_list.append(y_norm_i)
    y = y_norm_list
    rnge_y = [0,1.2]
    scale = 'linear' 
else:    
    y = np.abs(y_list)
    scale = 'log' 
    rnge_y = [np.min(y)/10,np.max(y)*10]
    
color = 'black'
title = 'Photon Spectrum'
x_title = '$\omega$'

y_title = 'S'
   
name = 'N=%d kd=%.2f $\Omega$=%.2f $\Delta$=%.2f'%(N,k0*d,Omega,Delta)
rnge_x = [np.min(omega_list),np.max(omega_list)]

file_name_curve = 'Photon_Spectrum_N=%d_kd=%.2f_Omega=%.2f_Delta=%.2f.pdf'%(N,k0*d,Omega,Delta)
file_name = path_plot + file_name_curve

plot.Plot(x,y,color,title,x_title,y_title,scale,N,k0,name,rnge_x,rnge_y,file_name)                

                

