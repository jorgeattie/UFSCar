from qutip import *
import numpy as np
import h5py 
import os 


import delta_function
import gamma_function
import hamiltonian
import lindbladian
import correlation_functions
import filename_generator

import position_vector
import Two_Dim_Maps

#General Parameters

Gamma = 1
N = 2
S = 2
F = N+S
level = 2
n_hat = [0,1,0]


#Geometry

k0=2*np.pi
k_vec = [0,k0,0]
k = np.linalg.norm(k_vec)
k_unity = np.divide(k_vec,k)
n_hat = np.divide(n_hat,np.linalg.norm(n_hat))
k_vec = k0*k_unity
d = 1/k0
p = [0,0,1]
P = [p]*N
R = position_vector.position_vector(d,N)

#Parameters of the system

g = delta_function.delta_escalar(N,R,k0,Gamma)
f = gamma_function.gamma_escalar(N,R,k0,Gamma)
Omega = 0.2*Gamma
Delta = 5

#Parameters of the sensors

epsilon = 0.05


#Parameters of the plot 
omega_points = 256
omega_list = np.linspace(-25,25,omega_points)

#Parameters of the dynamics

Gamma_s = Gamma
H_a = hamiltonian.Hamiltonian_gen(N,R,g,Omega,Delta,k_vec,level,F)
L_a = lindbladian.Lindbladian_gen(N,f,F,level)
i = 0
j = 0

#Parameters of the hdf document

path_data = "../results/data/QT/Matrix_hdf12_g2_N=2_kd=1_Omega=2_Delta=5_1/"
if not os.path.exists(path_data): 
    os.makedirs(path_data)

filename = path_data+filename_generator.filename_gen_sensor("I",N,S,k,d,Omega,Delta) 


#Options

g2 = 1
g3 = 0
g2_map = 0
g3_map = 0
G2_map = 0
G3_map = 0
g3_linear = 0
g3_linear_plot = 0
opt_points = 0
conditional_map = 0



#Codes
    
if g2:
    G = np.zeros((omega_points,omega_points))
    g = np.zeros((omega_points,omega_points))
    for omega_1 in omega_list:  
        for omega_2 in omega_list:            
            omega_list_sensor = [omega_1,omega_2]        
            H_s = hamiltonian.Hamiltonian_Sensors(N,k,n_hat,R,level,F,omega_list_sensor,epsilon)
            H = H_a + H_s       
            L_s = lindbladian.Lindbladian_Sensors(N,F,level,Gamma_s)
            c_ops = L_a + L_s            
            atol = 10**(-23)            
            result = steadystate(H,c_ops,tol = atol)            
            g2_s,G2_s = correlation_functions.g2_sensor(N,F,level,result)            
            print(omega_1,omega_2)            
            g[i,j] = g2_s 
            G[i,j] = G2_s
            j += 1
        j = 0
        i += 1
        
    f = h5py.File(filename, 'w')
    g1 = f.create_group('1')
    dset_g2 = g1.create_dataset('g2',(omega_points,omega_points),dtype = 'float64')    
    dset_G2 = g1.create_dataset('G2',(omega_points,omega_points),dtype = 'float64')
    dset_g2[:,:] = g[:,:]
    dset_G2[:,:] = G[:,:]
    
if g3:
    G = np.zeros((omega_points,omega_points))
    g = np.zeros((omega_points,omega_points))
    for omega_1 in omega_list:  
        for omega_2 in omega_list:            
            omega_list_sensor = [omega_1,omega_2,-omega_1 - omega_2]        
            H_s = hamiltonian.Hamiltonian_Sensors(N,k,n_hat,R,level,F,omega_list_sensor,epsilon)
            H = H_a + H_s       
            L_s = lindbladian.Lindbladian_Sensors(N,F,level,Gamma_s)
            c_ops = L_a + L_s            
            atol = 10**(-23)            
            result = steadystate(H,c_ops,tol = atol)            
            g3_s,G3_s = correlation_functions.g3_sensor(N,F,level,result)            
            print(g3_s,G3_s,omega_1,omega_2,-omega_1 - omega_2)            
            g[i,j] = g3_s 
            G[i,j] = G3_s
            j += 1
        j = 0
        i += 1
        
    f = h5py.File(filename, 'w')
    g1 = f.create_group('1')
    dset_g3 = g1.create_dataset('g3',(omega_points,omega_points),dtype = 'float64')    
    dset_G3 = g1.create_dataset('G3',(omega_points,omega_points),dtype = 'float64')
    dset_g3[:,:] = g[:,:]
    dset_G3[:,:] = G[:,:]
    
if g2_map:
    f = h5py.File(filename,'r')
    g2_Matrix = f['/1/g2']
    
    M = g2_Matrix
    x = omega_list
    y = omega_list
    x_title = '$\omega_{2}$'
    y_title = '$\omega_{1}$'
    file_name = 'g2_map_N=%d_kd=%.2f_Omega=%.2f.pdf'%(N,k0*d,Omega)
    title = '$g^{(2)}(0)$'
    
    Two_Dim_Maps.TwoDim_MapLogScale(M,N,k0,d,x,y,x_title,y_title,file_name,title)
        
if g3_map:
    f = h5py.File(filename,'r')
    g3_Matrix = f['/1/g3']
    
    M = g3_Matrix
    x = omega_list
    y = omega_list
    x_title = '$\omega_{2}$'
    y_title = '$\omega_{1}$'
    file_name = 'g3_map_N=%d_kd=%.2f_Omega=%.2f.pdf'%(N,k0*d,Omega)
    title = '$g^{(3)}(0)$'
    
    Two_Dim_Maps.TwoDim_MapLogScale(M,N,k0,d,x,y,x_title,y_title,file_name,title)


if G2_map:
    f = h5py.File(filename,'r')
    G2_Matrix = f['1/G2']
    
    M = G2_Matrix
    x = omega_list
    y = omega_list
    x_title = '$\omega_{2}$'
    y_title = '$\omega_{1}$'
    file_name = 'Intensity_2_map_N=%d_kd=%.2f_Omega=%.2f.pdf'%(N,k0*d,Omega)
    title = 'Intensity'
    Two_Dim_Maps.TwoDim_MapLogScale(M,N,k0,d,x,y,x_title,y_title,file_name,title)


if G3_map:
    f = h5py.File(filename,'r')
    G3_Matrix = f['1/G3']
    
    M = G3_Matrix
    x = omega_list
    y = omega_list
    x_title = '$\omega_{2}$'
    y_title = '$\omega_{1}$'
    file_name = 'Intensity_3_map_N=%d_kd=%.2f_Omega=%.2f.pdf'%(N,k0*d,Omega)
    title = '$G^{(3)}(0)$'
    Two_Dim_Maps.TwoDim_MapLogScale(M,N,k0,d,x,y,x_title,y_title,file_name,title)
        

if conditional_map:
    
    f = h5py.File(filename,'r')
    g2_Matrix = f['/1/g2']
    G2_Matrix = f['1/G2']
    
    conditional_Matrix = np.zeros((omega_points,omega_points))
    
    for i in range(0,omega_points):
        for j in range(0,omega_points):
            if g2_Matrix[i,j] > 10:
                conditional_Matrix[i,j] = G2_Matrix[i,j]
            else:
                conditional_Matrix[i,j] = float('nan')                
    
    
    
    
    M = conditional_Matrix
    x = omega_list
    y = omega_list
    title = 'Conditional Map'
    x_title = '$\omega_{2}$'
    y_title = '$\omega_{1}$'
    file_name = 'Conditional_map_N=%d_kd=%.2f_Omega=%.2f.pdf'%(N,k0*d,Omega)
    Two_Dim_Maps.TwoDim_MapLogScale(M,N,k0,d,x,y,x_title,y_title,file_name,title)
    
             
             
             
             



