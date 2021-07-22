#Lindbladian

from qutip import *
import numpy as np
import matplotlib.pyplot as plt
import scipy.linalg as la

import sigma_minus
import sigma_plus

def Lindbladian(N,f):
    L = []
    L_parc = 0
    for i in range(0,N):
        SgmM_i = sigma_minus.Sigmam(N,2,i)
        for j in range(0,N):
            SgmP_j = sigma_plus.Sigmap(N,2,j)
            L_parc += f[i,j]*(2*(qutip.spre(SgmM_i)*qutip.spost(SgmP_j)) - (qutip.spre(SgmP_j*SgmM_i)+qutip.spost(SgmP_j*SgmM_i)))    
        L.append(0.5*L_parc)
    #print(L)
    return L

def Lindbladian_gen(N,f,F,level):
    L = []
    L_parc = 0
    for i in range(0,N):
        SgmM_i = sigma_minus.Sigmam_gen(level,F,i)
        for j in range(0,N):
            SgmP_j = sigma_plus.Sigmap_gen(level,F,i)
            L.append(0.5*f[i,j]*(2*(qutip.spre(SgmM_i)*qutip.spost(SgmP_j)) - (qutip.spre(SgmP_j*SgmM_i) + qutip.spost(SgmP_j*SgmM_i))))    
            #L.append(0.5*L_parc)
    #print(L)
    return L

def Lindbladian_Sensors(N,F,level,Gamma_s):
    L = []
    L_parc = 0
    for i in range(N,F):
        ksiP = sigma_plus.Sigmap_gen(level,F,i)
        ksiM = sigma_minus.Sigmam_gen(level,F,i)
        L.append(0.5*Gamma_s*(2*(qutip.spre(ksiM)*qutip.spost(ksiP)) - qutip.spre(ksiP*ksiM) - qutip.spost(ksiP*ksiM)))
        #L.append((Gamma_s/2)*L_parc)
    return L

        
        
        
        
        
        
        
        
        
        
        
        