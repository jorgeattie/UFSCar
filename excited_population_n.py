# -*- coding: utf-8 -*-
"""
Created on Thu Apr  8 10:02:54 2021

@author: Jorge
"""

from itertools import combinations
from qutip import *
import numpy as np
import matplotlib.pyplot as plt
import scipy.linalg as la


def excited_population_n(N,n):
    LISTA = [i for i in range(0,N)]
    index = list(combinations(LISTA,n))
    list_state = [basis(2,1)]*N
    projector_total = 0
    for item in index:
        for i in range (0,n):
            x = item[i]
            list_state[x] = basis(2,0)        
        psi_n = (tensor(list_state)).unit()
        projector_total += psi_n*psi_n.dag()
        for j in range(0,n):
            x = item[j]
            list_state[x] = basis(2,1)
    return projector_total

