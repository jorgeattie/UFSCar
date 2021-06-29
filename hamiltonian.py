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

def Hamiltonian_gen(N,R,g,Omega,Delta,k_vec,level,F):
    H1,H2,H3 = 0,0,0
    for i in range(0,N):
        SgmP_i = sigma_plus.Sigmap_gen(level,F,i)
        SgmM_i = sigma_minus.Sigmam_gen(level,F,i)
        H1 += -Delta*(SgmP_i*SgmM_i)        
        H2 +=  -0.5*Omega*(np.exp(1j*np.dot(k_vec,R[i]))*SgmP_i + np.exp(-1j*np.dot(k_vec,R[i]))*SgmM_i)
        for j in range(0,N):
            SgmM_j = sigma_minus.Sigmam_gen(level,F,i)
            H3 += g[i,j]*SgmP_i*SgmM_j
    H = H1 + H2 + H3    
    return H

def Hamiltonian_Sensors(N,k,n_hat,R,level,F,Omega,epsilon):
    E_p = electric_field.E_gen(N,level,F,'plus',n_hat,R,k)
    E_m = electric_field.E_gen(N,level,F,'minus',n_hat,R,k)
    for i in range(N,F):
        ksiP = sigma_plus.Sigmap_gen(level,F,i)
        ksiM = sigma_minus.Sigmam_gen(level,F,i)
        H_s_1 = Omega*ksiP*ksiM
        H_s_2 = E_m*ksiM + E_p*ksiP
    H_s = H_s_1 + epsilon*H_s_2
    return H_s




