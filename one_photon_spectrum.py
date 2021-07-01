# -*- coding: utf-8 -*-
"""
Created on Fri Jun 25 14:35:38 2021

@author: Usuario
"""

from qutip import *
from scipy.fft import fft

import numpy as np
import matplotlib.pyplot as plt
import scipy.linalg as la


import position_vector
import correlation_functions
import delta_function
import gamma_function
import lindbladian
import plot
import hamiltonian
import excited_population_n

#General Parameters

Gamma = 1
N = 2
S = 2
level = 2
n_hat = [0,1,0]

#Geometry

k0=2*np.pi
k_vec = [0,k0,0]
k = np.linalg.norm(k_vec)
d = 0.05/k0
p = [0.577,0.577,0.577]
P = [p]*N
R = position_vector.position_vector(d,N)
g = delta_function.g_matrix(Gamma,R,P,k0,N)
f = gamma_function.f_matrix(Gamma,R,P,k0,N)
Delta = 0
Omega = 30



F = N + S
Gamma_s = Gamma
epsilon = 10**(-4)

time_points = 100
time_reference = 10
taulist = np.linspace(0,time_reference,time_points)

omega_points = 100
omega_list = np.linspace(-60,60,omega_points)
count = 0

n1_list = []
n_exc_1 = excited_population_n.excited_population_n_gen(N,1,F)

for omega_sensor in omega_list:
    H_gen = hamiltonian.Hamiltonian_gen(N,R,g,Omega,Delta,k_vec,level,F)
    H_s = hamiltonian.Hamiltonian_Sensors(N,k,n_hat,R,level,F,omega_sensor,epsilon)
    L_gen = lindbladian.Lindbladian_gen(N,f,F,level)
    L_s = lindbladian.Lindbladian_Sensors(N,F,level,Gamma_s)
    H = H_gen + H_s
    c_ops = L_gen + L_s
    result = steadystate(H,c_ops)
    n1 = expect(n_exc_1,result)
    n1_list.append(n1)
    count += 1
    print(count)
    
x = omega_list
y = n1_list
color = 'red'
title = 'Population'
x_title = '$\omega_{s}$'
y_title = '$P_{1}$'
scale = 'log'
name = 'Population for 1 excitation'
rnge = [10**(-23),10**(6)]
plot.Plot(x,y,color,title,x_title,y_title,scale,N,k0,d,name,rnge)



  









