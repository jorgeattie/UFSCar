from qutip import *
import numpy as np
import h5py 
import os 

import filename_generator
import position_vector
import hamiltonian
import delta_function
import lindbladian
import gamma_function
import correlation_functions

#General Parameters

Gamma = 1
Gamma_s = Gamma
omega_points = 256
omega_list = np.linspace(-25,25,omega_points)
epsilon = 0.05
Delta_list = [0,5,10]
N_list = [1,2]
Omega_list = [0.2*Gamma,2*Gamma,10*Gamma]
distance_list = [0.04,0.1,1,2*np.pi]
count_N1 = 0
count_N2 = 0
S = 2
level = 2
n_hat = [0,1,0]
vmin_g = 10**(-3)
vmax_g = 10**(3)
atol = 10**(-23)


for Delta in Delta_list:
    for N in N_list:                                                
        if N == 1:            
            F = N + S
            R = [[0,0,0]]
            k0=2*np.pi
            k_vec = [0,k0,0]
            k = np.linalg.norm(k_vec)
            d = 0
            g = delta_function.delta_escalar(N,R,k0,Gamma)
            f = gamma_function.gamma_escalar(N,R,k0,Gamma)
            L_s = lindbladian.Lindbladian_Sensors(N,F,level,Gamma_s)
            L_a = lindbladian.Lindbladian_gen(N,f,F,level)            
            c_ops = L_s + L_a            
            for Omega in Omega_list:
                
                path_data = "../results/data/QT/Matrix_hdf12_Bell_Inequalities_2_N=%d_Omega=%.2f_Delta=%.2f/"%(N,Omega,Delta)
                if not os.path.exists(path_data): 
                     os.makedirs(path_data)
                     
                filename = path_data+filename_generator.filename_gen_sensor("I",N,S,k,d,Omega,Delta) 
                                               
                i,j = 0,0
                B_s = np.zeros((omega_points,omega_points))
                H_a = hamiltonian.Hamiltonian_gen(N,R,g,Omega,Delta,k_vec,level,F)
                
                for omega_1 in omega_list:
                    for omega_2 in omega_list:                          
                        
                        omega_list_sensor = [omega_1,omega_1,omega_1,omega_1]
                        H_s = hamiltonian.Hamiltonian_Sensors(N,k,n_hat,R,level,F,omega_list_sensor,epsilon)
                        H = H_a + H_s                               
                        result = steadystate(H,c_ops,tol = atol)
                        g2_s_1111,b_1111 = correlation_functions.g2_sensor(N,F,level,result) 
                        
                        omega_list_sensor = [omega_2,omega_2,omega_2,omega_2]
                        H_s = hamiltonian.Hamiltonian_Sensors(N,k,n_hat,R,level,F,omega_list_sensor,epsilon)
                        H = H_a + H_s  
                        result = steadystate(H,c_ops,tol = atol)
                        g2_s_2222,b_2222 = correlation_functions.g2_sensor(N,F,level,result)
                        
                        omega_list_sensor = [omega_1,omega_2,omega_2,omega_1]
                        H_s = hamiltonian.Hamiltonian_Sensors(N,k,n_hat,R,level,F,omega_list_sensor,epsilon)
                        H = H_a + H_s 
                        result = steadystate(H,c_ops,tol = atol)
                        g2_s_1221,b_1221 = correlation_functions.g2_sensor(N,F,level,result)
                        
                        omega_list_sensor = [omega_1,omega_1,omega_2,omega_2]
                        H_s = hamiltonian.Hamiltonian_Sensors(N,k,n_hat,R,level,F,omega_list_sensor,epsilon)
                        H = H_a + H_s 
                        result = steadystate(H,c_ops,tol = atol)
                        g2_s_1122,b_1122 = correlation_functions.g2_sensor(N,F,level,result)
                        
                        omega_list_sensor = [omega_2,omega_2,omega_1,omega_1]
                        H_s = hamiltonian.Hamiltonian_Sensors(N,k,n_hat,R,level,F,omega_list_sensor,epsilon)
                        H = H_a + H_s 
                        result = steadystate(H,c_ops,tol = atol)
                        g2_s_2211,b_2211 = correlation_functions.g2_sensor(N,F,level,result)
                        
                        B_s_value = np.sqrt(2)*(np.abs(np.divide((b_1111 + b_2222 - 4*b_1221 - b_1122 - b_2211),(b_1111+b_2222 + 2*b_1221))))
                        
                        if B_s_value > 2:
                            B_s[i,j] = True
                        else:
                            B_s[i,j] = False
                              
                        print(N,omega_1,omega_2)
                        
                        j += 1                    
                    j = 0
                    i += 1
                    
                    
                f = h5py.File(filename, 'w')
                g1 = f.create_group('1')                  
                dset_Bs = g1.create_dataset('Bs',(omega_points,omega_points),dtype = 'float64')                
                dset_Bs[:,:] = B_s[:,:]
                f.close()
                
        
        if N == 2:
            
            F = N + S
            
            for distance in distance_list:                    
                k0=2*np.pi
                k_vec = [0,k0,0]
                k = np.linalg.norm(k_vec)
                k_unity = np.divide(k_vec,k)
                n_hat = np.divide(n_hat,np.linalg.norm(n_hat))
                k_vec = k0*k_unity
                d = distance/k0
                p = [0,0,1]
                P = [p]*N
                R = position_vector.position_vector(d,N)
                g = delta_function.delta_escalar(N,R,k0,Gamma)
                f = gamma_function.gamma_escalar(N,R,k0,Gamma)
                L_s = lindbladian.Lindbladian_Sensors(N,F,level,Gamma_s)
                L_a = lindbladian.Lindbladian_gen(N,f,F,level)            
                c_ops = L_s + L_a           
                
                for Omega in Omega_list:
                    path_data = "../results/data/QT/Matrix_hdf12_Bell_Inequalities__2_N=%d_kd=%.2f_Omega=%.2f_Delta=%.2f/"%(N,k0*d,Omega,Delta)
                    if not os.path.exists(path_data): 
                        os.makedirs(path_data)
                    
                    filename = path_data+filename_generator.filename_gen_sensor("BI",N,S,k,d,Omega,Delta)
                    
                    i,j = 0,0
                    B_s = np.zeros((omega_points,omega_points))
                    H_a = hamiltonian.Hamiltonian_gen(N,R,g,Omega,Delta,k_vec,level,F)
                    
                    for omega_1 in omega_list:
                        for omega_2 in omega_list:                        
                            
                            omega_list_sensor = [omega_1,omega_1,omega_1,omega_1]
                            H_s = hamiltonian.Hamiltonian_Sensors(N,k,n_hat,R,level,F,omega_list_sensor,epsilon)
                            H = H_a + H_s                               
                            result = steadystate(H,c_ops,tol = atol)
                            g2_s_1111,b_1111 = correlation_functions.g2_sensor(N,F,level,result) 
                        
                            omega_list_sensor = [omega_2,omega_2,omega_2,omega_2]
                            H_s = hamiltonian.Hamiltonian_Sensors(N,k,n_hat,R,level,F,omega_list_sensor,epsilon)
                            H = H_a + H_s  
                            result = steadystate(H,c_ops,tol = atol)
                            g2_s_2222,b_2222 = correlation_functions.g2_sensor(N,F,level,result)
                        
                            omega_list_sensor = [omega_1,omega_2,omega_2,omega_1]
                            H_s = hamiltonian.Hamiltonian_Sensors(N,k,n_hat,R,level,F,omega_list_sensor,epsilon)
                            H = H_a + H_s 
                            result = steadystate(H,c_ops,tol = atol)
                            g2_s_1221,b_1221 = correlation_functions.g2_sensor(N,F,level,result)
                        
                            omega_list_sensor = [omega_1,omega_1,omega_2,omega_2]
                            H_s = hamiltonian.Hamiltonian_Sensors(N,k,n_hat,R,level,F,omega_list_sensor,epsilon)
                            H = H_a + H_s 
                            result = steadystate(H,c_ops,tol = atol)
                            g2_s_1122,b_1122 = correlation_functions.g2_sensor(N,F,level,result)
                        
                            omega_list_sensor = [omega_2,omega_2,omega_1,omega_1]
                            H_s = hamiltonian.Hamiltonian_Sensors(N,k,n_hat,R,level,F,omega_list_sensor,epsilon)
                            H = H_a + H_s 
                            result = steadystate(H,c_ops,tol = atol)
                            g2_s_2211,b_2211 = correlation_functions.g2_sensor(N,F,level,result)
                        
                            B_s_value = np.sqrt(2)*(np.abs(np.divide((b_1111 + b_2222 - 4*b_1221 - b_1122 - b_2211),(b_1111+b_2222 + 2*b_1221))))
                        
                            if B_s_value > 2:
                                B_s[i,j] = True
                            else:
                                B_s[i,j] = False
                        
                            j += 1
                        j = 0
                        i += 1
                        print(N,omega_1,omega_2)
                    
                    f = h5py.File(filename, 'w')
                    g1 = f.create_group('1')                  
                    dset_Bs = g1.create_dataset('Bs',(omega_points,omega_points),dtype = 'float64')                
                    dset_Bs[:,:] = B_s[:,:]
                    f.close()
                    
                                       
                    
                    
                    
                    
            
            
            
                
                        
                        
                        
                                        
                                        
                                        
                                        

                
                                                                                                    