#Sigma plus

from qutip import *
import numpy as np
import matplotlib.pyplot as plt
import scipy.linalg as la

def Sigmap(N,level,i):
    q_list = [qeye(level)]*N
    q_list[i] = sigmap()
    sigmapN = tensor(q_list)
    return sigmapN

def Sigmap_gen(level,F,i):
    q_list = [qeye(level)]*F
    q_list[i] = sigmap()
    sigmapN = tensor(q_list)
    return sigmapN

