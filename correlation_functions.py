# -*- coding: utf-8 -*-
"""
Created on Thu Apr  8 15:25:59 2021

@author: Jorge
"""

from qutip import *
import numpy as np
import matplotlib.pyplot as plt
import scipy.linalg as la

import hamiltonian
import electric_field
import sigma_plus
import sigma_minus



def g2_0(N,R,g,Omega,Delta,k_vec,c_ops,k,n_hat,level):
    H = hamiltonian.Hamiltonian(N,R,g,Omega,Delta,k_vec)
    result = steadystate(H,c_ops)
    E_p = electric_field.E(N,k,n_hat,R,level,'plus')
    E_m = electric_field.E(N,k,n_hat,R,level,'minus')
    numerator = expect(E_m*E_m*E_p*E_p,result)
    denominator = (expect(E_m*E_p,result))**2
    g2_0 = np.divide(numerator,denominator)
    G2_0 = numerator
    return g2_0,G2_0

def g3_0(N,R,g,Omega,Delta,k_vec,c_ops,k,n_hat,level):
    H = hamiltonian.Hamiltonian(N,R,g,Omega,Delta,k_vec)
    result = steadystate(H,c_ops)
    E_p = electric_field.E(N,k,n_hat,R,level,'plus')
    E_m = electric_field.E(N,k,n_hat,R,level,'minus')
    numerator = expect(E_m*E_m*E_m*E_p*E_p*E_p,result)
    denominator = (expect(E_m*E_p,result))**3
    g3_0 = np.divide(numerator,denominator)
    return g3_0

def g4_0(N,R,g,Omega,Delta,k_vec,c_ops,k,n_hat,level):
    H = hamiltonian.Hamiltonian(N,R,g,Omega,Delta,k_vec)
    result = steadystate(H,c_ops)
    E_p = electric_field.E(N,k,n_hat,R,level,'plus')
    E_m = electric_field.E(N,k,n_hat,R,level,'minus')
    numerator = expect(E_m*E_m*E_m*E_m*E_p*E_p*E_p*E_p,result)
    denominator = (expect(E_m*E_p,result))**4
    g4_0 = np.divide(numerator,denominator)
    return g4_0


def G1_0(N,R,g,Omega,Delta,k_vec,c_ops,k,n_hat,level):
    H = hamiltonian.Hamiltonian(N,R,g,Omega,Delta,k_vec)
    result = steadystate(H,c_ops)
    E_p = electric_field.E(N,k,n_hat,R,level,'plus')
    E_m = electric_field.E(N,k,n_hat,R,level,'minus')
    G1_0 = expect(E_m*E_p,result)
    return G1_0


def G3_0(N,R,g,Omega,Delta,k_vec,c_ops,k,n_hat,level):
    H = hamiltonian.Hamiltonian(N,R,g,Omega,Delta,k_vec)
    result = steadystate(H,c_ops)
    E_p = electric_field.E(N,k,n_hat,R,level,'plus')
    E_m = electric_field.E(N,k,n_hat,R,level,'minus')
    G3_0 = expect(E_m*E_m*E_m*E_p*E_p*E_p,result)    
    return G3_0

def G4_0(N,R,g,Omega,Delta,k_vec,c_ops,k,n_hat,level):
    H = hamiltonian.Hamiltonian(N,R,g,Omega,Delta,k_vec)
    result = steadystate(H,c_ops)
    E_p = electric_field.E(N,k,n_hat,R,level,'plus')
    E_m = electric_field.E(N,k,n_hat,R,level,'minus')
    G4_0 = expect(E_m*E_m*E_m*E_m*E_p*E_p*E_p*E_p,result)    
    return G4_0

def g1_T(H,c_ops,N,k,n_hat,R,level,taulist,F):
    result = steadystate(H,c_ops)
    E_p = electric_field.E_gen(N,level,F,'plus',n_hat,R,k)
    E_m = electric_field.E_gen(N,level,F,'minus',n_hat,R,k)
    num = correlation_ss(H,taulist,c_ops,E_p,E_m,solver = 'me',reverse=True)
    den = expect(E_m*E_p,result)
    g1_T = num/den
    return g1_T

def g2_T(N,R,g,Omega,Delta,k_vec,c_ops,k,n_hat,level,taulist):
     H = hamiltonian.Hamiltonian(N,R,g,Omega,Delta,k_vec)
     result = steadystate(H,c_ops)
     E_p = electric_field.E(N,k,n_hat,R,level,'plus')
     E_m = electric_field.E(N,k,n_hat,R,level,'minus')
     G2_T = correlation_4op_1t(H,result,taulist,c_ops,E_m,E_m,E_p,E_p,solver = 'me')
     I2_T = correlation_ss(H,taulist,c_ops,E_m*E_p,1,solver='me')
     I2 = expect(E_m*E_p,result)    
     g2_T = np.divide(G2_T,I2*I2_T)
     return g2_T,G2_T,I2_T,I2


def g3_T(N,R,g,Omega,Delta,k_vec,c_ops,k,n_hat,level,taulist):
     H = hamiltonian.Hamiltonian(N,R,g,Omega,Delta,k_vec)
     result = steadystate(H,c_ops)
     E_p = electric_field.E(N,k,n_hat,R,level,'plus')
     E_m = electric_field.E(N,k,n_hat,R,level,'minus')
     G3_T = correlation_3op_1t(H,result,taulist,c_ops,E_m*E_m,E_m*E_p,E_p*E_p,solver = 'me')
     I3 = expect(E_m*E_p,result)     
     I3_T = correlation_ss(H,taulist,c_ops,E_m*E_p,1,solver='me')
     g3_T = G3_T/((I3**2)*I3_T)
     return g3_T,G3_T,I3_T,I3

def g21_T(N,R,g,Omega,Delta,k_vec,c_ops,k,n_hat,level,taulist):
     H = hamiltonian.Hamiltonian(N,R,g,Omega,Delta,k_vec)
     result = steadystate(H,c_ops)
     E_p = electric_field.E(N,k,n_hat,R,level,'plus')
     E_m = electric_field.E(N,k,n_hat,R,level,'minus')
     A = E_m*E_m
     B = E_m*E_p
     C = E_p*E_p
     g21_T = np.divide(correlation_3op_1t(H,result,taulist,c_ops,A,B,C,solver = 'me'),
                       np.multiply(expect(E_m*E_m*E_p*E_p,result),expect(E_m*E_p,result)))
     return g21_T
 

def intensity_sensor(N,F,level,result):
    ksiP_list = []
    ksiM_list = []
    for i in range(N,F):
        ksiP_i = sigma_plus.Sigmap_gen(level,F,i)
        ksiM_i = sigma_minus.Sigmam_gen(level,F,i)
        ksiP_list.append(ksiP_i)
        ksiM_list.append(ksiM_i)
    intensity = expect(ksiP_list[0]*ksiP_list[1]*ksiM_list[0]*ksiM_list[1],result)
    return intensity
    

def g2_sensor(N,F,level,result):  
    ksiP_list = []
    ksiM_list = []
    for i in range(N,F):
        ksiP_i = sigma_plus.Sigmap_gen(level,F,i)
        ksiM_i = sigma_minus.Sigmam_gen(level,F,i)
        ksiP_list.append(ksiP_i)
        ksiM_list.append(ksiM_i)
    num = expect(ksiP_list[0]*ksiP_list[1]*ksiM_list[0]*ksiM_list[1],result)
    den = expect(ksiP_list[0]*ksiM_list[0],result)*expect(ksiP_list[1]*ksiM_list[1],result)
    g2 = num/den
    G2 = num
    return g2,G2

def g3_sensor(N,F,level,result):  
    ksiP_list = []
    ksiM_list = []
    for i in range(N,F):
        ksiP_i = sigma_plus.Sigmap_gen(level,F,i)
        ksiM_i = sigma_minus.Sigmam_gen(level,F,i)
        ksiP_list.append(ksiP_i)
        ksiM_list.append(ksiM_i)
    num = expect(ksiP_list[0]*ksiP_list[1]*ksiP_list[2]*ksiM_list[0]*ksiM_list[1]*ksiM_list[2],result)
    den = expect(ksiP_list[0]*ksiM_list[0],result)*expect(ksiP_list[1]*ksiP_list[2]*ksiM_list[1]*ksiM_list[2],result)
    g3 = num/den
    G3 = num
    return g3,G3


    
        
        
        
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    


