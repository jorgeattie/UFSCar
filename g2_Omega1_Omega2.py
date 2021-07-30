from qutip import *
import numpy as np
import matplotlib.pyplot as plt
import scipy.linalg as la
import h5py 
import os 
import matplotlib as mpl
from matplotlib.cm import datad, get_cmap
import matplotlib.colors as mcolors
from matplotlib import rc

import delta_function
import gamma_function
import hamiltonian
import lindbladian
import correlation_functions
import filename_generator


def rgb_color(a,b,c):
    norm = np.float(255)
    a = a/norm
    b = b/norm
    c = c/norm
    return (a,b,c)

def make_colormap(seq):
    """Return a LinearSegmentedColormap
    seq: a sequence of floats and RGB-tuples. The floats should be increasing
    and in the interval (0,1).
    """
    seq = [(None,) * 3, 0.0] + list(seq) + [1.0, (None,) * 3]
    cdict = {'red': [], 'green': [], 'blue': []}
    for i, item in enumerate(seq):
        if isinstance(item, float):
            r1, g1, b1 = seq[i - 1]
            r2, g2, b2 = seq[i + 1]
            cdict['red'].append([item, r1, r2])
            cdict['green'].append([item, g1, g2])
            cdict['blue'].append([item, b1, b2])
    return mcolors.LinearSegmentedColormap('CustomMap', cdict)

#General Parameters

Gamma = 1
N = 2
S = 2
level = 2
n_hat = [0.3,0.1,-0.5]


#Geometry

k0=2*np.pi
k_vec = [0,k0,0]
k = np.linalg.norm(k_vec)
k_unity = np.divide(k_vec,k)
n_hat = np.divide(n_hat,np.linalg.norm(n_hat))
k_vec = k0*k_unity
d = (0.05*2*np.pi)/k0
#d = 0.05/k0
p = [0,0,1]
P = [p]*N
R = [[0,0,0],np.multiply(d,[1,0,0])]
#R = position_vector.position_vector(d,N)
g = delta_function.g_matrix(Gamma,R,P,k0,N)
f = gamma_function.f_matrix(Gamma,R,P,k0,N)
Omega = 30*Gamma
Delta = 0
epsipower = 4
epsilon = 10**(-epsipower)  

run = 0
plot = 1

path_data = "../results/data/QT/Matrix_hdf5_7/" 
if not os.path.exists(path_data): 
    os.makedirs(path_data)

filename = path_data+filename_generator.filename_gen("g2",N,S,k,d,Omega,Delta) 



if run: 
    #Sensor Parameters
    F = N+S    
    omega_points = 512
    omega_list = np.linspace(-70,70,omega_points)
    g2_s_Omega_1_list = []
    g2_s_Omega_2_list = []
    Gamma_s = Gamma
    i = 0
    j = 0
    S_i = 0
    g2_s_Matrix = np.zeros((omega_points,omega_points))
    H_a = hamiltonian.Hamiltonian_gen(N,R,g,Omega,Delta,k_vec,level,F)
    L_a = lindbladian.Lindbladian_gen(N,f,F,level)
    
    
    for omega_1 in omega_list:  
        for omega_2 in omega_list: 
            print(omega_1,omega_2)        
            omega_list_sensor = [omega_1,omega_2]        
            H_s = hamiltonian.Hamiltonian_Sensors(N,k,n_hat,R,level,F,omega_list_sensor,epsilon)
            H = H_a + H_s       
            L_s = lindbladian.Lindbladian_Sensors(N,F,level,Gamma_s)
            c_ops = L_a + L_s
            result = steadystate(H,c_ops)
            g2_s = correlation_functions.g2_sensor(N,F,level,result)
            g2_s_Matrix[i,j]  = g2_s
            j += 1
        j = 0
        i += 1
            
    #dset_t = g1.create_dataset("t", (len(t_list),), dtype='float64')
   
    f = h5py.File(filename, 'w')
    g1 = f.create_group('1')
    
    dset_g2 = g1.create_dataset('g2',(omega_points,omega_points),dtype = 'float64')
    dset_g2[:,:] = g2_s_Matrix[:,:]
    dset_omegaList = g1.create_dataset("omegaList",(len(omega_list),),dtype = 'float64')
    dset_omegaList[:] = np.array(omega_list)[:]
    
    
if plot:    

    f = h5py.File(filename,'r')
    g2_s_Matrix = f['/1/g2']
    omega_list = f['/1/omegaList']


    c1, c2, c3 = rgb_color(73, 143, 201), rgb_color(255,255,255), rgb_color(204, 95, 73)
    custom_cmap = make_colormap([c1,c2,0.15,c2,0.16,c2,0.17,c2,c3])

    extent = np.min(omega_list), np.max(omega_list), np.min(omega_list), np.max(omega_list)
    print(np.min(g2_s_Matrix), np.max(g2_s_Matrix))
    plt.imshow(np.log10(np.abs(g2_s_Matrix)),cmap=custom_cmap,extent = extent, origin="lower", vmin=np.log10(0.4),vmax=np.log10(35),interpolation="spline36") #,vmin=np.log10(0.4),vmax=np.log10(35))
    #fig, ax = plt.subplots(figsize=(1,1))
    #norm = mcolors.LogNorm(vmin = np.min(g2_s_Matrix),vmax = np.max(g2_s_Matrix))
    norm = mpl.colors.Normalize(vmin = np.min(g2_s_Matrix),vmax = np.max(g2_s_Matrix))
    #cb = mpl.colorbar.ColorbarBase(ax,cmap=custom_cmap,norm=norm,orientation='vertical')
    plt.colorbar(mpl.cm.ScalarMappable(norm=norm,cmap=custom_cmap))
    plt.savefig("g2_w1_w2_kd=%.2f_Delta=%.1f_epspower=%d.pdf" %(k0*d,Delta,epsipower) ,bbox_inches="tight")
    #plt.show()    
        
        
        
    



