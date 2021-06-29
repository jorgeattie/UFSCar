#Delta Function
from qutip import *
import numpy as np
import matplotlib.pyplot as plt
import scipy.linalg as la

import green_tensor

def g_matrix(Gamma,R,P,k0,N):
    g = np.add(np.zeros((N,N)),np.multiply(1j,np.zeros((N,N))))
    for i in range(0,N):
        for j in range(0,N):
            G_ij = green_tensor.Green_Tensor(Gamma,R,k0,i,j)
            G_real = np.real(G_ij)
            g[i,j] = np.tensordot(np.conj(P[i]),np.tensordot(G_real,P[j],axes=1),axes=1)   
    return g