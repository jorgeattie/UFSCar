# -*- coding: utf-8 -*-
"""
Created on Fri Mar  5 16:01:14 2021

@author: Jorge
"""


from qutip import *
import numpy as np
import matplotlib.pyplot as plt
import scipy.linalg as la


def diag_function(M):
    eigs = la.eig(M)
    eigvls = eigs[0]
    eigvls_order = sorted(eigvls)
    Gamma_diag = np.real(eigvls_order[-1]) 
    return Gamma_diag

