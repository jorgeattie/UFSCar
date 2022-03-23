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

def Sigmam_gen_mLevel(N,S,F,level_a,level_s,i):
    q_list_1 = [qeye(level_a)]*N
    q_list_2 = [qeye(level_s)]*S
    q_list = q_list_1 + q_list_2
    
    if i >= 0 and i < N:
        q_list[i] = destroy(level_a)
    if i >= N and i < F:
        q_list[i] = destroy(level_s)
        
    sigmamN = tensor(q_list)
    return sigmamN

def Sigmam_atoms(N,level_a,i):
    q_list = [qeye(level_a)]*N
    q_list[i] = sigmam()
    sigmamN_atoms = tensor(q_list)
    return sigmamN_atoms
    
def Sigmam_sensors(S,level_s,i):
    q_list = [qeye(level_s)]*S
    q_list[i] = destroy(level_s)
    sigmamN_sensors = tensor(q_list)
    return sigmamN_sensors
    

    


