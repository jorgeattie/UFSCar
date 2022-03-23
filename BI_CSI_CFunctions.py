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
import sigma_plus
import sigma_minus

#General Parameters

Gamma = 1
Gamma_s = Gamma
omega_points = 128
omega_list = np.linspace(-25,25,omega_points)
epsilon = 0.05
Delta_list = [0*Gamma]
N_list = [1]
Omega_list = [0.1*Gamma,1*Gamma,5*Gamma]
#distance_list = [0.04,0.1,1,2*np.pi]
count_N1 = 0
count_N2 = 0
S = 2
level_a = 2
level_s = 3
n_hat = [0,1,0]
atol = 10**(-23)

B_s_value_list = []
R_cs_list = []

B_s = np.zeros((omega_points,omega_points))
R_cs = np.zeros((omega_points,omega_points))
g2_Matrix = np.zeros((omega_points,omega_points))

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
            L_s = lindbladian.Lindbladian_Sensors(N,F,level_a,level_s,Gamma_s,S)
            L_a = lindbladian.Lindbladian_gen(N,f,F,level_a,level_s,S)            
            c_ops = L_s + L_a            
            for Omega in Omega_list:
                
                #path_data = "../results/data/QT/BI_FF/Matrix_hdf12_Bell_Inequalities_Fixed_Frequencies_N=%d_S=%d_Omega=%.2f_Delta=%.2f/"%(N,S,Omega,Delta)
                #if not os.path.exists(path_data): 
                #     os.makedirs(path_data)
                     
                #filename = path_data+filename_generator.filename_gen_sensor("BI",N,S,k0,d,Omega,Delta) 
                                               
                i,j = 0,0                
                B_s = np.zeros((omega_points,omega_points))
                H_a = hamiltonian.Hamiltonian_gen(N,R,g,Omega,Delta,k_vec,level_a,level_s,F,S)                
                
                for omega_1 in omega_list:
                    for omega_2 in omega_list:
                        omega_sensor_list = [omega_1,omega_2]                
                        H_s = hamiltonian.Hamiltonian_Sensor_BI_FF(N,level_a,level_s,n_hat,R,k,F,omega_sensor_list,epsilon,S)              
                        H = H_a + H_s                            
                        atol = 10**(-23)
                        result = steadystate(H,c_ops,tol = atol)                                                                                                                         
                
                        ksiP_list = []
                        ksiM_list = []
                        for k in range(N,F):
                            ksiP_i = sigma_plus.Sigmap_gen_mLevel(N,S,F,level_a,level_s,i)
                            ksiM_i = sigma_minus.Sigmam_gen_mLevel(N,S,F,level_a,level_s,i)                          
                            ksiP_list.append(ksiP_i)
                            ksiM_list.append(ksiM_i)
                        
                        b_1 = expect(ksiP_list[0]*ksiP_list[0]*ksiM_list[0]*ksiM_list[0],result)
                        b_2 = expect(ksiP_list[1]*ksiP_list[1]*ksiM_list[1]*ksiM_list[1],result)
                        b_3 = expect(ksiP_list[0]*ksiP_list[1]*ksiM_list[1]*ksiM_list[0],result)
                        b_4 = expect(ksiP_list[0]*ksiP_list[0]*ksiM_list[1]*ksiM_list[1],result)
                        b_5 = expect(ksiP_list[1]*ksiP_list[1]*ksiM_list[0]*ksiM_list[0],result)
               
                        
                        num = b_1 + b_2 -4*b_3 - b_4 - b_5
                        den = b_1 + b_2 + 2*b_3
                        
                        # Correlation Functions
                        
                        g_12 = (b_3)/((expect(ksiP_list[0]*ksiM_list[0],result))*(expect(ksiP_list[1]*ksiM_list[1],result)))
                        g_11 = (b_1)/((expect(ksiP_list[0]*ksiM_list[0],result))*(expect(ksiP_list[0]*ksiM_list[0],result)))
                        g_22 = (b_2)/((expect(ksiP_list[1]*ksiM_list[1],result))*(expect(ksiP_list[1]*ksiM_list[1],result)))
                        
                        g2_Matrix[i,j] = g_12
                        
                        
                        #Cauchy-Schwarz inequality
                        
                        R_cs_value = (g_12**2)/(g_11*g_22)
                        R_cs_list.append(R_cs_value)
                        R_cs[i,j] = R_cs_value        

                                                                                                                          
                
                        Bell inequality
                        
                        B_s_value = np.sqrt(2)*np.abs(num/den)
                        B_s_value_list.append(B_s_value)
                        B_s[i,j] = B_s_value
                        
                        print()                                                                                                                                                                        
                        
                        j += 1                    
                    j = 0
                    i += 1
                    
                    print(omega_1,omega_2)                                                           
                
        
      
                                        
                                        
                                        

                
                                                                                                    


