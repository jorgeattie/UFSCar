#Sigma Minus

from qutip import *
import numpy as np
import matplotlib.pyplot as plt
import scipy.linalg as la

def Sigmam(N,level,i):
    q_list = [qeye(level)]*N
    q_list[i] = sigmam()
    sigmamN = tensor(q_list)
    return sigmamN

def Sigmam_gen(level,F,i):
    q_list = [qeye(level)]*F
    q_list[i] = sigmam()
    sigmamN = tensor(q_list)
    return sigmamN

#F = N + S

    


