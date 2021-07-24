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

#General Parameters

Gamma = 1
N = 2
S = 1
level = 2
n_hat = [0,1,0]

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
Delta_s = Gamma*g[0,1]
Delta_a = -Gamma*g[0,1]

Gamma_s = Gamma
epsilon = 10**(-4)

omega_points = 500
omega_list = np.linspace(-80,80,omega_points)
n_list = []
F = N+S
index_sensor = 0


ksiM = sigma_minus.Sigmam_gen(level,F,N+index_sensor)
ksiP = sigma_plus.Sigmap_gen(level,F,N+index_sensor)
n_sensor = ksiP*ksiM  


for omega_sensor in omega_list:
     print(omega_sensor)     
     H_gen = hamiltonian.Hamiltonian_gen(N,R,g,Omega,Delta_a,k_vec,level,F)
     H_s = hamiltonian.Hamiltonian_Sensors(N,k,n_hat,R,level,F,omega_sensor,epsilon)
     H = H_gen + H_s
     L_gen = lindbladian.Lindbladian_gen(N,f,F,level)
     L_s = lindbladian.Lindbladian_Sensors(N,F,level,Gamma_s)
     c_ops = L_gen + L_s
     result = steadystate(H,c_ops)            
     n_sensor_i = expect(n_sensor,result)
     n_list.append(n_sensor_i)
                      
x = omega_list
y = n_list
color = 'red'
title = '1PS - One photon spectrum'
x_title = '$\omega$'
y_title = '$S$'
scale = 'log'
name = "$N = %d \quad kd = %.2f \quad \Omega = %.1f \quad \Delta = %.1f$" %(N,k0*d,Omega,Delta_a)
rnge = [0.8*np.min(n_list),1.1*np.max(n_list)]

plot.Plot(x,y,color,title,x_title,y_title,scale,N,k0,d,name,rnge)




























#time_points = 1000
#time_reference = 10
#taulist = np.linspace(0,time_reference,time_points)

#H_gen = hamiltonian.Hamiltonian_gen(N,R,g,Omega,Delta,k_vec,level,F)
#H_s = hamiltonian.Hamiltonian_Sensors(N,k,n_hat,R,level,F,omega_sensor,epsilon)
#H = H_gen + H_s
#L_gen = lindbladian.Lindbladian_gen(N,f,F,level)
#L_s = lindbladian.Lindbladian_Sensors(N,F,level,Gamma_s)
#c_ops = L_gen + L_s
#g1 = correlation_functions.g1_T(H,c_ops,N,k,n_hat,R,level,taulist,F)
#S = fft(g1)
#amp = np.abs(S)

#freq_points = 1000
#freq_list = np.linspace(-60,60,freq_points)


#x = freq_list
#y = amp
#color = 'red'
#title = '1PS - One photon spectrum'
#x_title = '$\omega$'
#y_title = '$S$'
#scale = 'log'
#name = "$N = %d \quad kd = %.2f \quad \Omega = %.1f \quad \Delta = %.1f$" %(N,k0*d,Omega,Delta)
#rnge = [10**(-4),10**(4)]

#plot.Plot(x,y,color,title,x_title,y_title,scale,N,k0,d,name,rnge)





    
    




