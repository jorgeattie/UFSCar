# -*- coding: utf-8 -*-
"""
Created on Fri Jun 25 14:35:38 2021

@author: Usuario
"""

from qutip import *
from scipy.fft import fft
from scipy.fft import fftfreq

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
Omega = 30
omega_sensor = 100
Delta = -1200



F = N + S
Gamma_s = Gamma
epsilon = 10**(-4)

time_points = 1000
time_reference = 10
taulist = np.linspace(0,time_reference,time_points)

H_gen = hamiltonian.Hamiltonian_gen(N,R,g,Omega,Delta,k_vec,level,F)
H_s = hamiltonian.Hamiltonian_Sensors(N,k,n_hat,R,level,F,omega_sensor,epsilon)
H = H_gen + H_s
L_gen = lindbladian.Lindbladian_gen(N,f,F,level)
L_s = lindbladian.Lindbladian_Sensors(N,F,level,Gamma_s)
c_ops = L_gen + L_s
g1 = correlation_functions.g1_T(H,c_ops,N,k,n_hat,R,level,taulist,F)
S = fft(g1)
amp = np.abs(S)

freq_points = 1000
freq_list = np.linspace(-60,60,freq_points)


x = freq_list
y = amp
color = 'red'
title = '1PS - One photon spectrum'
x_title = '$\omega$'
y_title = '$S$'
scale = 'log'
name = "$N = %d \quad kd = %.2f \quad \Omega = %.1f \quad \Delta = %.1f$" %(N,k0*d,Omega,Delta)
rnge = [10**(-4),10**(4)]

plot.Plot(x,y,color,title,x_title,y_title,scale,N,k0,d,name,rnge)





    
    




