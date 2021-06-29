#Gamma Function

from qutip import *
import numpy as np
import matplotlib.pyplot as plt
import scipy.linalg as la

import green_tensor

def f_matrix(Gamma,R,P,k0,N):
    f = np.add(np.zeros((N,N)),np.multiply(1j,np.zeros((N,N))))
    for i in range(0,N):
        for j in range(0,N):
            G_ij = green_tensor.Green_Tensor(Gamma,R,k0,i,j)
            G_Img = np.imag(G_ij)
            f[i,j]=np.tensordot(np.conj(P[i]),np.tensordot(np.multiply(2,G_Img),P[j],axes=1),axes=1)
    return f