from qutip import *
import numpy as np
import h5py 
import os 


import hamiltonian
import lindbladian
import position_vector
import delta_function
import gamma_function
import filename_generator
import sigma_plus
import sigma_minus
import electric_field


#General Parameters

Gamma = 1
Gamma_s = Gamma
omega_points = 1024
omega_list = np.linspace(-200,200,omega_points)
epsilon = 10**(-2)
Delta = 0
N = 2
S = 2
Omega = 30*Gamma
level = 2
n_hat = [0,1,0]
atol = 10**(-23)
k0=2*np.pi
k_vec = [0,k0,0]
k = np.linalg.norm(k_vec)
d = 0.05/k
i,j = 0,0

scalar_model = 1
vector_model = 0
sensor_method = 1
spectrum = 1

ksiP_list = []
ksiM_list = []
P_list = []
P_norm_list = []
H_s,L_s = 0,[]


#Parameters of the atomic system
R = position_vector.position_vector(d,N)
p = [0,0,1]
P = [p]*N

#Models

if vector_model:
    g = delta_function.g_matrix(Gamma,R,P,k0,N)
    f = gamma_function.f_matrix(Gamma,R,P,k0,N)

if scalar_model:        
    g = delta_function.delta_escalar(N,R,k0,Gamma)
    f = gamma_function.gamma_escalar(N,R,k0,Gamma)
    
#Parameter for the populations
F = N + S

#hdf files
if sensor_method:
    path_data = "../results/data/QT/Photon_Spectrum/Matrix_PSCT_WS_hdf5_N=%d_S=%d_Omega=%.2f_Delta=%.2f_epsilon=%.2f/"%(N,S,Omega,Delta,epsilon)   
    if not os.path.exists(path_data): 
        os.makedirs(path_data)
else:
    path_data = "../results/data/QT/Photon_Spectrum/Matrix_PSCT_hdf5_N=%d_S=%d_Omega=%.2f_Delta=%.2f_epsilon=%.2f/"%(N,S,Omega,Delta,epsilon)   
    if not os.path.exists(path_data): 
        os.makedirs(path_data)            

filename = path_data+filename_generator.filename_gen_sensor("PSCT",N,S,k0,d,Omega,Delta)     

     

  
if sensor_method:     
    H_a = hamiltonian.Hamiltonian_gen(N,R,g,Omega,Delta,k_vec,level,F,S)
    L_a = lindbladian.Lindbladian_gen(N,f,F,level,S)
    L_s = lindbladian.Lindbladian_Sensors(N,F,level,Gamma_s,S)
    L = L_a + L_s
    if spectrum:
        H_s = 0    
        for omega in omega_list:
            print(omega)
            omega_sensor_list = [omega]*S
            H_s += hamiltonian.Hamiltonian_Sensors(N,k,n_hat,R,level,F,omega_sensor_list,epsilon)        
        H = H_a + H_s                      
        E_m = electric_field.E_gen(N,level,F,'minus',n_hat,R,k)
        E_p = electric_field.E_gen(N,level,F,'plus',n_hat,R,k)
        S = spectrum_ss(H,omega_list,L,E_m,E_p)
        y_list = S       
    else:
        for i in range(N,F):
            ksiP_i = sigma_plus.Sigmap_gen(level,F,i)
            ksiM_i = sigma_minus.Sigmam_gen(level,F,i)
            ksiP_list.append(ksiP_i)
            ksiM_list.append(ksiM_i)
        for omega in omega_list:
            print(omega)
            omega_sensor_list = [omega]*S
            H_s = hamiltonian.Hamiltonian_Sensors(N,k,n_hat,R,level,F,omega_sensor_list,epsilon)
            H = H_a + H_s
            result = steadystate(H,L,tol = 10**(-23))
            P_i = expect(ksiP_list[0]*ksiM_list[0]*ksiP_list[1]*ksiM_list[1],result)
            P_list.append(P_i)
        y_list = P_list                                                                    
else:        
    H = hamiltonian.Hamiltonian(N,R,g,Omega,Delta,k_vec)
    L = lindbladian.Lindbladian(N,f)
    E_m = electric_field.E(N,k,n_hat,R,level,'minus')
    E_p = electric_field.E(N,k,n_hat,R,level,'plus')
    S = spectrum_ss(H,omega_list,L,E_m,E_p)
    y_list = S

f = h5py.File(filename, 'w')
g1 = f.create_group('1')                  
dset_y_list_data = g1.create_dataset('y_list_data',((omega_points),),dtype = 'float64')                
dset_y_list_data[:] = y_list[:]
f.close()
    









    

