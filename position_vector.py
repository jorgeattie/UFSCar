#Position vector

from qutip import *
import numpy as np
import matplotlib.pyplot as plt
import scipy.linalg as la

def position_vector(d,N):
    R = [[0,0,0]]
    j = 0
    while j < (N-1):
        j = j + 1
        R.append([j*d,0,0])
    print(R)
    return R

def position_vector_two_two(dist,d):
    if dist == 'far':
        R = [[0,0,0],[d,0,0],[100,0,0],[100 + d,0,0]]        
    else:
        R = [[0,0,0],[d,0,0],[d+2*d,0,0],[d+2*d+d,0,0]]
        
    print(R)
    return R

def position_vector_two_one(dist,d):
    if dist == 'far':
        R = [[0,0,0],[d,0,0],[100 + (d/2),0,0]]
    else:
        R = [[0,0,0],[d,0,0],[d+2*d,0,0]]
    
    print(R)
    return R

    
