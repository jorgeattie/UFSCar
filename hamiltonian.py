#Hamiltonian

from qutip import *
import numpy as np
import matplotlib.pyplot as plt
import scipy.linalg as la

import sigma_minus
import sigma_plus
import electric_field

def Hamiltonian(N,R,g,Omega,Delta,k_vec):
    H1,H2,H3 = 0,0,0
    for i in range(0,N):
        SgmP_i = sigma_plus.Sigmap(N,2,i)
        SgmM_i = sigma_minus.Sigmam(N,2,i)
        H1 += -Delta*(SgmP_i*SgmM_i)
        
        H2 +=  -0.5*Omega*(np.exp(1j*np.dot(k_vec,R[i]))*SgmP_i + np.exp(-1j*np.dot(k_vec,R[i]))*SgmM_i)
        for j in range(0,N):
            SgmM_j = sigma_minus.Sigmam(N,2,j)
            H3 += g[i,j]*SgmP_i*SgmM_j
    H = H1 + H2 + H3    
    return H

def Hamiltonian_gen(N,R,g,Omega,Delta,k_vec,level_a,level_s,F,S):
    H1,H2,H3 = 0,0,0
    for i in range(0,N):
        SgmP_i = sigma_plus.Sigmap_gen_mLevel(N,S,F,level_a,level_s,i)
        SgmM_i = sigma_minus.Sigmam_gen_mLevel(N,S,F,level_a,level_s,i)
        H1 += -Delta*(SgmP_i*SgmM_i)        
        H2 +=  -0.5*Omega*(np.exp(1j*np.dot(k_vec,R[i]))*SgmP_i + np.exp(-1j*np.dot(k_vec,R[i]))*SgmM_i)
        for j in range(0,N):
            SgmM_j = sigma_minus.Sigmam_gen_mLevel(N,S,F,level_a,level_s,i)
            H3 += g[i,j]*SgmP_i*SgmM_j
    H = H1 + H2 + H3    
    return H

def Hamiltonian_Sensors(N,k,n_hat,R,level,F,omega_sensor_list,epsilon):
    E_p = electric_field.E_gen(N,level,F,'plus',n_hat,R,k)
    E_m = electric_field.E_gen(N,level,F,'minus',n_hat,R,k)
    H_s_1 = 0
    H_s_2 = 0
    index_omega = 0
    for i in range(N,F):
        omega_sensor = omega_sensor_list[index_omega]
        index_omega +=1
        ksiP = sigma_plus.Sigmap_gen(level,F,i)
        ksiM = sigma_minus.Sigmam_gen(level,F,i)
        H_s_1 += omega_sensor*ksiP*ksiM
        H_s_2 += E_m*ksiM + E_p*ksiP
    H_s = H_s_1 + epsilon*H_s_2
    return H_s

def Hamiltonian_Sensor_Frequency_Filtered_Correlation_Two_Sensors(N,level,F,n_hat,R,k,omega_sensor_list,b1111_label,b2222_label,b1221_label,b1122_label,b2211_label,epsilon):
    E_p = electric_field.E_gen(N,level,F,'plus',n_hat,R,k)
    E_m = electric_field.E_gen(N,level,F,'minus',n_hat,R,k)    
    ksiP_list = []
    ksiM_list = []
    
    for i in range(N,F):        
        ksiP = sigma_plus.Sigmap_gen(level,F,i)
        ksiM = sigma_minus.Sigmam_gen(level,F,i)
        ksiP_list.append(ksiP)
        ksiM_list.append(ksiM)
    
    if b1111_label == 1:
        ksiP_list[1] = ksiP_list[0]
        ksiM_list[1] = ksiM_list[0]
        H_s = omega_sensor_list[0]*ksiP_list[0]*ksiM_list[0] + omega_sensor_list[1]*ksiP_list[1]*ksiM_list[1]
        + epsilon*(E_m*ksiM_list[0] + E_p*ksiP_list[0] + E_m*ksiM_list[1] + E_p*ksiP_list[1])
    
    if b2222_label == 1:
        ksiP_list[0] = ksiP_list[1]
        ksiM_list[0] = ksiM_list[1]
        H_s = omega_sensor_list[0]*ksiP_list[0]*ksiM_list[0] + omega_sensor_list[1]*ksiP_list[1]*ksiM_list[1]
        + epsilon*(E_m*ksiM_list[0] + E_p*ksiP_list[0] + E_m*ksiM_list[1] + E_p*ksiP_list[1])
        
    if b1221_label == 1:
        H_s = omega_sensor_list[0]*ksiP_list[0]*ksiM_list[0] + omega_sensor_list[1]*ksiP_list[1]*ksiM_list[1]
        + epsilon*(E_m*ksiM_list[0] + E_p*ksiP_list[0] + E_m*ksiM_list[1] + E_p*ksiP_list[1])
    
    if b1122_label == 1:
        ksiP_list[1] = ksiP_list[0]
        ksiM_list[0] = ksiM_list[1]
        H_s = omega_sensor_list[0]*ksiP_list[0]*ksiM_list[0] + omega_sensor_list[1]*ksiP_list[1]*ksiM_list[1]
        + epsilon*(E_m*ksiM_list[0] + E_p*ksiP_list[0] + E_m*ksiM_list[1] + E_p*ksiP_list[1])
    
    if b2211_label == 1:
        ksiP_list[0] = ksiP_list[1]
        ksiM_list[1] = ksiM_list[0]
        H_s = omega_sensor_list[0]*ksiP_list[0]*ksiM_list[0] + omega_sensor_list[1]*ksiP_list[1]*ksiM_list[1]
        + epsilon*(E_m*ksiM_list[0] + E_p*ksiP_list[0] + E_m*ksiM_list[1] + E_p*ksiP_list[1])
    
    return H_s

def Hamiltonian_Sensor_BI_FF(N,level_a,level_s,n_hat,R,k,F,omega_sensor_list,epsilon,S):
    E_p = electric_field.E_gen_mLevel(N,'plus',n_hat,R,k,S,F,level_a,level_s)
    E_m = electric_field.E_gen_mLevel(N,'minus',n_hat,R,k,S,F,level_a,level_s)
    H_s_1 = 0
    H_s_2 = 0
    index_omega = 0    
    for i in range(N,F):
        omega_sensor = omega_sensor_list[index_omega]
        index_omega +=1
        ksiP = sigma_plus.Sigmap_gen_mLevel(N,S,F,level_a,level_s,i)
        ksiM = sigma_minus.Sigmam_gen_mLevel(N,S,F,level_a,level_s,i)
        #print(E_m,E_p,ksiM,ksiP)
        #E_m_conv,E_p_conv,ksiM_conv,ksiP_conv = np.array(E_m),np.array(E_p),np.array(ksiM),np.array(ksiP)
        H_s_1 += omega_sensor*ksiP*ksiM          
        H_s_2 += E_m*ksiM + E_p*ksiP                   
        #H_s_1 = omega_sensor*ksiP_conv*ksiM_conv
        #H_s_2 ++ E_m_conv*ksiM_conv + E_p_conv*ksiP_conv
    H_s = H_s_1 + epsilon*H_s_2
    return H_s


def Hamiltonian_Sensor_BI(N,level,F,n_hat,R,k,omega_sensor_list,epsilon):
     E_p = electric_field.E_gen(N,level,F,'plus',n_hat,R,k)
     E_m = electric_field.E_gen(N,level,F,'minus',n_hat,R,k)     
     ksiP_1 = sigma_plus.Sigmap_gen(level,F,1) 
     ksiP_2 = sigma_plus.Sigmap_gen(level,F,2)
     ksiM_1 = sigma_minus.Sigmam_gen(level,F,1)
     ksiM_2 = sigma_minus.Sigmam_gen(level,F,2)
     H_s_1 = omega_sensor_list[0]*ksiP_1*ksiM_1 + omega_sensor_list[1]*ksiP_2*ksiM_2
     H_s_2 = E_m*ksiM_1 + E_p*ksiP_1 + E_m*ksiM_2 + E_p*ksiP_2
     H_s = H_s_1 + epsilon*H_s_2
     return H_s
         
         
         
     
    




