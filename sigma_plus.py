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

def Sigmap_gen_mLevel(N,S,F,level_a,level_s,i):
    q_list_1 = [qeye(level_a)]*N
    q_list_2 = [qeye(level_s)]*S
    q_list = q_list_1 + q_list_2
    
    if i >= 0 and i < N:
        q_list[i] = sigmap()
    if i >= N and i < F:
        q_list[i] = create(level_s)
        
    sigmapN = tensor(q_list)
    return sigmapN


def Sigmap_atoms(N,level_a,i):
    q_list = [qeye(level_a)]*N
    q_list[i] = sigmap()
    sigmapN_atoms = tensor(q_list)
    return sigmapN_atoms
    
def Sigmap_sensors(S,level_s,i):
    q_list = [qeye(level_s)]*S
    q_list[i] = create(level_s)
    sigmapN_sensors = tensor(q_list)
    return sigmapN_sensors
    
