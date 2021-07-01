# -*- coding: utf-8 -*-
"""
Created on Thu May  6 14:07:36 2021

@author: Jorge
"""

from qutip import *
import numpy as np
import matplotlib.pyplot as plt
import scipy.linalg as la
from itertools import combinations

import h5py 
import os 
import filename_generator

import position_vector
import Two_Dim_Maps

path_data_QT = "../results/data/QT/Matrix_hdf5_11/" 

#General Parameters

Gamma = 1
N = 4
level = 2
n_hat = [0,1,0]

#Geometry

k0=2*np.pi
k_vec = [0,k0,0]
k = np.linalg.norm(k_vec)
d = 0.1/k0
p = [0,0,1]
P = [p]*N
R = position_vector.position_vector(d,N)

#Plot parameters
len_Delta = 100
len_Omega = 100
points_delta = np.linspace(-4000,4000,len_Delta)
points_omega = np.linspace(10**(-10),2000,len_Omega)


filename_QT = path_data_QT + filename_generator.filename_gen('QT',N,k0,d)
f = h5py.File(filename_QT,'r')
P1_matrix = f['/1/P1']
P2_matrix = f['/1/P2']
G2_matrix = f['/1/g2']
G3_matrix = f['/1/g3']
G4_matrix = f['/1/g4']

M = P2_matrix[:,:]/((P1_matrix[:,:])**2)
x = points_omega
y = points_delta
x_title = '$\Omega$'
y_title = '$\Delta$'

Two_Dim_Maps.TwoDim_Map(M,N,k0,d,x,y,x_title,y_title)














