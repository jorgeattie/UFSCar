#Green Tensor

from qutip import *
import numpy as np
import matplotlib.pyplot as plt
import scipy.linalg as la

def Green_Tensor(Gamma,R,k,i,j):
    G = np.add(np.zeros((3,3)),np.multiply(1j,np.zeros((3,3))))
    if (i != j):        
        R_ij = np.add(R[i],np.multiply(-1.0,R[j]))
        r_ij = np.linalg.norm(R_ij)
        G = (3.0*Gamma/4.0)*((np.exp(1j*k*r_ij))/(((k*r_ij)**3)))*((k*r_ij)**2 + 1j*k*r_ij -1)*np.identity(3)-((k*r_ij)**2 + 1j*3*k*r_ij -3)*np.tensordot(R_ij,R_ij,axes=0)
                                     
    else:
        G = 1j*Gamma/2*np.identity(3)
    return G