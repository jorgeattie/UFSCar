# -*- coding: utf-8 -*-
"""
Created on Fri Mar  5 15:38:11 2021

@author: Jorge
"""

from qutip import *
import numpy as np
import matplotlib.pyplot as plt
import scipy.linalg as la

import sigma_plus
import sigma_minus

def E(N,k,n_hat,R,level,signal):
    E = 0
    for i in range(0,N):
        if signal == 'plus':
            E += np.exp(-1j*k*np.dot(n_hat,R[i]))*sigma_minus.Sigmam(N,level,i)
        else:
            E +=  np.exp(1j*k*np.dot(n_hat,R[i]))*sigma_plus.Sigmap(N,level,i)   
    return E

def E_gen(N,level,F,signal,n_hat,R,k):
    E = 0
    for i in range(0,N):
        if signal == 'plus':
            E += np.exp(-1j*k*np.dot(n_hat,R[i]))*sigma_minus.Sigmam_gen(level,F,i)
        else:
            E +=  np.exp(1j*k*np.dot(n_hat,R[i]))*sigma_plus.Sigmap_gen(level,F,i)   
    return E

