from qutip import *
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
import sigma_minus
import sigma_plus 
import electric_field

#General Parameters

Gamma = 1
N = 2
S = 1
level = 2
n_hat = [0,0,1]

#Geometry

k0=2*np.pi
k_vec = [0,k0,0]
k = np.linalg.norm(k_vec)
k_unity = np.divide(k_vec,k)
n_hat = np.divide(n_hat,np.linalg.norm(n_hat))
k_vec = k0*k_unity
d = 0.5/k0
p = [0,0,1]
P = [p]*N
R = [[0,0,0],np.multiply(d/np.sqrt(3),[1,1,1])]
#R = position_vector.position_vector(d,N)
g = delta_function.g_matrix(Gamma,R,P,k0,N)
f = gamma_function.f_matrix(Gamma,R,P,k0,N)
Omega = 30*Gamma
#Delta_s = Gamma*g[0,1]
#Delta_a = -Gamma*g[0,1]
Delta_a = 0

H = hamiltonian.Hamiltonian(N,R,g,Omega,Delta_a,k_vec)
c_ops = lindbladian.Lindbladian(N,f)
omega_points = 500
wlist = np.linspace(-80,80,omega_points)
E_p = electric_field.E(N,k,n_hat,R,level,'plus')
E_m = electric_field.E(N,k,n_hat,R,level,'minus')

S = spectrum_ss(H,wlist,c_ops,E_m,E_p)

x = wlist
y = S
color = 'red'
title = '1PS - One photon spectrum'
x_title = '$\omega$'
y_title = '$S$'
scale = 'log'
name = "$N = %d \quad kd = %.2f \quad \Omega = %.1f \quad \Delta = %.1f$" %(N,k0*d,Omega,Delta_a)
rnge = [0.8*np.min(S),1.1*np.max(S)]

plot.Plot(x,y,color,title,x_title,y_title,scale,N,k0,d,name,rnge)
















