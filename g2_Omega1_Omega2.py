from qutip import *
import numpy as np
import matplotlib.pyplot as plt
import scipy.linalg as la
import h5py 
import os 


import delta_function
import gamma_function
import hamiltonian
import lindbladian
import correlation_functions
import filename_generator

#General Parameters

Gamma = 1
N = 2
S = 2
level = 2
n_hat = [0,-1,1]

#Geometry

k0=2*np.pi
k_vec = [0,k0,0]
k = np.linalg.norm(k_vec)
k_unity = np.divide(k_vec,k)
n_hat = np.divide(n_hat,np.linalg.norm(n_hat))
k_vec = k0*k_unity
d = (0.05)/k0
#d = 0.05/k0
p = [0,0,1]
P = [p]*N
R = [[0,0,0],np.multiply(d/np.sqrt(3),[1,1,1])]
#R = position_vector.position_vector(d,N)
g = delta_function.delta_escalar(N,R,k0,Gamma)
f = gamma_function.gamma_escalar(N,R,k0,Gamma)
Omega = 30*Gamma
Delta = 0

#Sensor Parameters
F = N+S
epsilon = 10**(-4)
omega_points = 256
omega_list = np.linspace(-70,70,omega_points)
g2_s_Omega_1_list = []
g2_s_Omega_2_list = []
Gamma_s = Gamma
i = 0
j = 0
S_i = 0
g2_s_Matrix = np.zeros((omega_points,omega_points))
H_a = hamiltonian.Hamiltonian_gen(N,R,g,Omega,Delta,k_vec,level,F)
L_a = lindbladian.Lindbladian_gen(N,f,F,level)

for omega_1 in omega_list:  
    for omega_2 in omega_list: 
        print(omega_1,omega_2)        
        omega_list_sensor = [omega_1,omega_2]        
        H_s = hamiltonian.Hamiltonian_Sensors(N,k,n_hat,R,level,F,omega_list_sensor,epsilon)
        H = H_a + H_s       
        L_s = lindbladian.Lindbladian_Sensors(N,F,level,Gamma_s)
        c_ops = L_a + L_s
        result = steadystate(H,c_ops)
        g2_s = correlation_functions.g2_sensor(N,F,level,result)
        g2_s_Matrix[i,j]  = g2_s
        j += 1
    j = 0
    i += 1
    
#path_data_QT = "../results/data/QT/Matrix_hdf5/" 
#if not os.path.exists(path_data_QT): 
#    os.makedirs(path_data_QT)

#filename_QT = path_data_QT+filename_generator.filename_gen("QT",N,k0,d)
#f = h5py.File(filename_QT, 'w')
#g1 = f.create_group('1')

#dset_g2 = g1.create_dataset('g2',(omega_points,omega_points),dtype = 'float64')
#dset_g2[:,:] = g2_s_Matrix[:,:]

extent = np.min(omega_list), np.max(omega_list), np.min(omega_list), np.max(omega_list)
plt.imshow(np.log10(np.abs(g2_s_Matrix)),cmap='coolwarm',extent = extent) #,vmin=np.log(0.4),vmax=np.log(30))
plt.colorbar()
    

        
        
        
    



