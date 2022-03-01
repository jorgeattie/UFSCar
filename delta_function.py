#Delta Function
from qutip import *
import numpy as np
import matplotlib.pyplot as plt
import scipy.linalg as la
import math

import green_tensor

def g_matrix(Gamma,R,P,k0,N):
    g = np.add(np.zeros((N,N)),np.multiply(1j,np.zeros((N,N))))
    for i in range(0,N):
        for j in range(0,N):
            G_ij = green_tensor.Green_Tensor(Gamma,R,k0,i,j)
            G_real = np.real(G_ij)
            g[i,j] = np.tensordot(np.conj(P[i]),np.tensordot(G_real,P[j],axes=1),axes=1)   
    return g

def delta_escalar(N,R,k0,Gamma):
    delta_matrix = np.zeros((N,N))
    for i in range(0,N):
        for j in range(0,N):
            if i != j:
                R_ij = np.add(R[i],np.multiply(-1.0,R[j]))
                r_ij = np.linalg.norm(R_ij)
                delta_matrix[i,j] = -0.5*Gamma*math.cos(k0*r_ij)/(k0*r_ij)
    return delta_matrix
                
            