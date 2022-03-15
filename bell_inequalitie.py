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
import sigma_plus
import sigma_minus

#General Parameters

Gamma = 1
Gamma_s = Gamma
omega_points = 10
omega_list = np.linspace(-25,25,omega_points)
epsilon = 20
Delta_list = [0,5,10,30]
N_list = [1,2]
Omega_list = [0.2*Gamma,2*Gamma,10*Gamma]
distance_list = [0.04,0.1,1,2*np.pi]
count_N1 = 0
count_N2 = 0
S = 2
level = 2
n_hat = [0,1,0]
atol = 10**(-23)


#Teste

num_N1_list = []
num_N2_list = []
den_N1_list = []
den_N2_list = []


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
                
                path_data = "../results/data/QT/Matrix_hdf12_Bell_Inequalities_N=%d_S=%d_Omega=%.2f_Delta=%.2f/"%(N,S,Omega,Delta)
                if not os.path.exists(path_data): 
                     os.makedirs(path_data)
                     
                filename = path_data+filename_generator.filename_gen_sensor("BI",N,S,k0,d,Omega,Delta) 
                                               
                i,j = 0,0                
                B_s = np.zeros((omega_points,omega_points))
                H_a = hamiltonian.Hamiltonian_gen(N,R,g,Omega,Delta,k_vec,level,F)                
                
                for omega_1 in omega_list:
                    for omega_2 in omega_list:
                        #omega_sensor_list = [omega_1,omega_2]                
                        #H_s = hamiltonian.Hamiltonian_Sensors(N,k,n_hat,R,level,F,omega_sensor_list,epsilon)                
                        #H = H_a + H_s                            
                        #atol = 10**(-23)
                        #result = steadystate(H,c_ops,tol = atol)  
                        
                        omega_sensor_list = [omega_1,omega_1]
                        H_s = hamiltonian.Hamiltonian_Sensor_BI(N,level,F,n_hat,R,k,omega_sensor_list,epsilon)
                        H = H_a + H_s
                        atol = 10**(-23)
                        result = steadystate(H,c_ops,tol = atol)                                                                        
                
                        ksiP_list = []
                        ksiM_list = []
                        for k in range(N,F):
                            ksiP_i = sigma_plus.Sigmap_gen(level,F,k)
                            ksiM_i = sigma_minus.Sigmam_gen(level,F,k)                            
                            ksiP_list.append(ksiP_i)
                            ksiM_list.append(ksiM_i)
                        
                        b11 = expect(ksiP_list[0]*ksiP_list[1]*ksiM_list[1]*ksiM_list[0],result)
                        
                        p_11_1 = expect(ksiP_list[0]*ksiM_list[0],result)
                        p_11_2 = expect(ksiP_list[1]*ksiM_list[1],result)
                        print(N,omega_1,omega_2)
                        print("Population S1: ",p_11_1)
                        print("Population S2: ",p_11_2)
                        
                        omega_sensor_list = [omega_2,omega_2]
                        H_s = hamiltonian.Hamiltonian_Sensor_BI(N,level,F,n_hat,R,k,omega_sensor_list,epsilon)
                        H = H_a + H_s
                        atol = 10**(-23)
                        result = steadystate(H,c_ops,tol = atol)                                                                        
                
                        ksiP_list = []
                        ksiM_list = []
                        for k in range(N,F):
                            ksiP_i = sigma_plus.Sigmap_gen(level,F,k)
                            ksiM_i = sigma_minus.Sigmam_gen(level,F,k)                            
                            ksiP_list.append(ksiP_i)
                            ksiM_list.append(ksiM_i)
                        
                        b22 = expect(ksiP_list[0]*ksiP_list[1]*ksiM_list[1]*ksiM_list[0],result)
                        
                        omega_sensor_list = [omega_1,omega_2]
                        H_s = hamiltonian.Hamiltonian_Sensor_BI(N,level,F,n_hat,R,k,omega_sensor_list,epsilon)
                        H = H_a + H_s
                        atol = 10**(-23)
                        result = steadystate(H,c_ops,tol = atol)                                                                        
                
                        ksiP_list = []
                        ksiM_list = []
                        for k in range(N,F):
                            ksiP_i = sigma_plus.Sigmap_gen(level,F,k)
                            ksiM_i = sigma_minus.Sigmam_gen(level,F,k)                            
                            ksiP_list.append(ksiP_i)
                            ksiM_list.append(ksiM_i)
                        
                        b12_1 = expect(ksiP_list[0]*ksiP_list[1]*ksiM_list[1]*ksiM_list[0],result)
                        b12_2 = expect(ksiP_list[0]*ksiP_list[0]*ksiM_list[1]*ksiM_list[1],result)
                                                
                        omega_sensor_list = [omega_2,omega_1]                        
                        H_s = hamiltonian.Hamiltonian_Sensor_BI(N,level,F,n_hat,R,k,omega_sensor_list,epsilon)
                        H = H_a + H_s
                        atol = 10**(-23)
                        result = steadystate(H,c_ops,tol = atol)                                                                        
                
                        ksiP_list = []
                        ksiM_list = []
                        for k in range(N,F):
                            ksiP_i = sigma_plus.Sigmap_gen(level,F,k)
                            ksiM_i = sigma_minus.Sigmam_gen(level,F,k)                            
                            ksiP_list.append(ksiP_i)
                            ksiM_list.append(ksiM_i)
                        
                        b21 = expect(ksiP_list[1]*ksiP_list[1]*ksiM_list[0]*ksiM_list[0],result)
                        
                        num = b11 + b22 -4*b12_1 - b12_2 - b21
                        den = b11 + b22 + 2*b12_1
                       
                        
                
                        #num = (expect(((ksiP_list[0]**2)*(ksiM_list[0]**2)),result) + expect(((ksiP_list[1]**2)*(ksiM_list[1]**2)),result) 
                        #- 4*expect((ksiP_list[0]*ksiP_list[1]*ksiM_list[1]*ksiM_list[0]),result) - expect((ksiP_list[0]**2)*(ksiM_list[1]**2),result) 
                        #- expect(((ksiP_list[1]**2)*(ksiM_list[0]**2)),result))
                        
                        
                        
                        #den = (expect(((ksiP_list[0]**2)*(ksiM_list[0]**2)),result) +  expect(((ksiP_list[1]**2)*(ksiM_list[1]**2)),result) 
                        #+ 2*expect((ksiP_list[0]*ksiP_list[1]*ksiM_list[1]*ksiM_list[0]),result))  
                        
                       

                        #num_N1_list.append(num)                    
                        #den_N1_list.append(den)
                                
                
                        B_s_value = np.sqrt(2)*np.abs(num/den)
                        #print(B_s_value,num/den)
                                                                             
                        
                        if B_s_value > 2:
                            B_s[i,j] = True
                        else:
                            B_s[i,j] = False
                              
                        #print(N,omega_1,omega_2)
                        
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
                
                    path_data = "../results/data/QT/Matrix_hdf12_Bell_Inequalities_N=%d_S=%d_Omega=%.2f_Delta=%.2f_kd=%.2f/"%(N,S,Omega,Delta,k0*d)
                    if not os.path.exists(path_data): 
                        os.makedirs(path_data)
                     
                    filename = path_data+filename_generator.filename_gen_sensor("BI",N,S,k0,d,Omega,Delta) 
                                               
                    i,j = 0,0                
                    B_s = np.zeros((omega_points,omega_points))
                    H_a = hamiltonian.Hamiltonian_gen(N,R,g,Omega,Delta,k_vec,level,F)                
                
                    for omega_1 in omega_list:
                        for omega_2 in omega_list:
                            #omega_sensor_list = [omega_1,omega_2]                
                            #H_s = hamiltonian.Hamiltonian_Sensors(N,k,n_hat,R,level,F,omega_sensor_list,epsilon)                
                            #H = H_a + H_s                            
                            #atol = 10**(-23)
                            #result = steadystate(H,c_ops,tol = atol)  
                        
                            omega_sensor_list = [omega_1,omega_1]
                            H_s = hamiltonian.Hamiltonian_Sensor_BI(N,level,F,n_hat,R,k,omega_sensor_list,epsilon)
                            H = H_a + H_s
                            atol = 10**(-23)
                            result = steadystate(H,c_ops,tol = atol)                                                                        
                            
                            ksiP_list = []
                            ksiM_list = []
                            for k in range(N,F):
                                ksiP_i = sigma_plus.Sigmap_gen(level,F,k)
                                ksiM_i = sigma_minus.Sigmam_gen(level,F,k)                            
                                ksiP_list.append(ksiP_i)
                                ksiM_list.append(ksiM_i)
                        
                            b11 = expect(ksiP_list[0]*ksiP_list[1]*ksiM_list[1]*ksiM_list[0],result)
                            p_11_1 = expect(ksiP_list[0]*ksiM_list[0],result)
                            p_11_2 = expect(ksiP_list[1]*ksiM_list[1],result)
                            print(N,omega_1,omega_2)
                            print("Population S1: ",p_11_1)
                            print("Population S2: ",p_11_2)
                        
                            omega_sensor_list = [omega_2,omega_2]
                            H_s = hamiltonian.Hamiltonian_Sensor_BI(N,level,F,n_hat,R,k,omega_sensor_list,epsilon)
                            H = H_a + H_s
                            atol = 10**(-23)
                            result = steadystate(H,c_ops,tol = atol)                                                                        
                            
                            ksiP_list = []
                            ksiM_list = []
                            for k in range(N,F):
                                ksiP_i = sigma_plus.Sigmap_gen(level,F,k)
                                ksiM_i = sigma_minus.Sigmam_gen(level,F,k)                            
                                ksiP_list.append(ksiP_i)
                                ksiM_list.append(ksiM_i)
                        
                            b22 = expect(ksiP_list[0]*ksiP_list[1]*ksiM_list[1]*ksiM_list[0],result)
                        
                            omega_sensor_list = [omega_1,omega_2]
                            H_s = hamiltonian.Hamiltonian_Sensor_BI(N,level,F,n_hat,R,k,omega_sensor_list,epsilon)
                            H = H_a + H_s
                            atol = 10**(-23)
                            result = steadystate(H,c_ops,tol = atol)                                                                        
                            
                            ksiP_list = []
                            ksiM_list = []
                            for k in range(N,F):
                                ksiP_i = sigma_plus.Sigmap_gen(level,F,k)
                                ksiM_i = sigma_minus.Sigmam_gen(level,F,k)                            
                                ksiP_list.append(ksiP_i)
                                ksiM_list.append(ksiM_i)
                        
                            b12_1 = expect(ksiP_list[0]*ksiP_list[1]*ksiM_list[1]*ksiM_list[0],result)
                            b12_2 = expect(ksiP_list[0]*ksiP_list[0]*ksiM_list[1]*ksiM_list[1],result)
                                                
                            omega_sensor_list = [omega_2,omega_1]                        
                            H_s = hamiltonian.Hamiltonian_Sensor_BI(N,level,F,n_hat,R,k,omega_sensor_list,epsilon)
                            H = H_a + H_s
                            atol = 10**(-23)
                            result = steadystate(H,c_ops,tol = atol)                                                                        
                
                            ksiP_list = []
                            ksiM_list = []
                            for k in range(N,F):
                                ksiP_i = sigma_plus.Sigmap_gen(level,F,k)
                                ksiM_i = sigma_minus.Sigmam_gen(level,F,k)                            
                                ksiP_list.append(ksiP_i)
                                ksiM_list.append(ksiM_i)
                        
                            b21 = expect(ksiP_list[1]*ksiP_list[1]*ksiM_list[0]*ksiM_list[0],result)
                        
                            num = b11 + b22 -4*b12_1 - b12_2 - b21
                            den = b11 + b22 + 2*b12_1
                       
                        
                
                            #num = (expect(((ksiP_list[0]**2)*(ksiM_list[0]**2)),result) + expect(((ksiP_list[1]**2)*(ksiM_list[1]**2)),result) 
                            #- 4*expect((ksiP_list[0]*ksiP_list[1]*ksiM_list[1]*ksiM_list[0]),result) - expect((ksiP_list[0]**2)*(ksiM_list[1]**2),result) 
                            #- expect(((ksiP_list[1]**2)*(ksiM_list[0]**2)),result))
                        
                        
                        
                            #den = (expect(((ksiP_list[0]**2)*(ksiM_list[0]**2)),result) +  expect(((ksiP_list[1]**2)*(ksiM_list[1]**2)),result) 
                            #+ 2*expect((ksiP_list[0]*ksiP_list[1]*ksiM_list[1]*ksiM_list[0]),result))  
                        
                       

                            #num_N1_list.append(num)                    
                            #den_N1_list.append(den)
                                
                
                            B_s_value = np.sqrt(2)*np.abs(num/den)
                            #print(B_s_value,num/den)
                                                                             
                        
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
                
                
                
                
               
                        
                           
                                
                    
                    
                    
                    
            
            
            
                
                        
                        
                        
                                        
                                        
                                        
                                        

                
                                                                                                    