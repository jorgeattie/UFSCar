from qutip import *
import numpy as np
import h5py 
import os 

import filename_generator
import position_vector
import Two_Dim_Maps

#General Parameters

Gamma = 1
Gamma_s = Gamma
omega_points = 256
omega_list = np.linspace(-25,25,omega_points)
epsilon = 0.05
Delta_list = [0,5,10]
N_list = [1,2]
Omega_list = [0.2*Gamma,2*Gamma,10*Gamma]
distance_list = [0.04,0.1,1,2*np.pi]
count_N1 = 0
count_N2 = 0
S = 2
level = 2
n_hat = [0,1,0]
vmin_g = 10**(-3)
vmax_g = 10**(3)


for Delta in Delta_list:
    for N in N_list:                                
        path_plot = "../plots/data/maps/Map_Correlation_Intensity_Conditional_CSI_N=%d"%(N)
        if not os.path.exists(path_plot): 
                    os.makedirs(path_plot)
        if N == 1:            
            F = N + S
            R = [[0,0,0]]
            k0=2*np.pi
            k_vec = [0,k0,0]
            k = np.linalg.norm(k_vec)
            d = 0
           
            for Omega in Omega_list:
                
                path_data = "../results/data/QT/Matrix_hdf12_Correlation_Intensity_Conditional_N=%d_Omega=%.2f_Delta=%.2f/"%(N,Omega,Delta)
                
                filename = path_data+filename_generator.filename_gen_sensor("I",N,S,k,d,Omega,Delta)
                f = h5py.File(filename,'r')
                g2_Matrix = f['/1/g2']

                R = np.zeros((omega_points,omega_points))
                conditional_Matrix = np.zeros((omega_points,omega_points))
                count = 0
                count_list = []
    
                for i in range(0,omega_points):
                    for j in range(0,omega_points):
                        R[i,j] = ((g2_Matrix[i,j])**(2))/(g2_Matrix[i,i]*g2_Matrix[j,j])

                g2_Condition = np.array(list(map(lambda x: x>10,g2_Matrix)))
                R_Condition = np.array(list(map(lambda y: y>1,R)))                
                Conditional_Matrix = np.zeros((omega_points,omega_points))
                M = np.zeros((omega_points,omega_points))

                for i in range(0,omega_points):
                    for j in range(0,omega_points):
                        Conditional_Matrix[i,j] = g2_Condition[i,j] and R_Condition[i,j]
                        if Conditional_Matrix[i,j] == True:
                            M[i,j] = g2_Matrix[i,j]
                        else: 
                            M[i,j] = float('nan')
                
                vmin,vmax = vmin_g,vmax_g
                x = omega_list
                y = omega_list
                title = 'Conditional Map - CSI'
                x_title = '$\omega_{2}$'
                y_title = '$\omega_{1}$'
                file_name_map = 'Conditional_map_CSI_N=%d_kd=%.2f_Omega=%.2f_Delta=%.2f.pdf'%(N,k0*d,Omega,Delta)
                file_name = path_plot + file_name_map
                Two_Dim_Maps.TwoDim_MapLogScale(M,N,k0,d,x,y,x_title,y_title,file_name,title,vmin,vmax)
                count_N1 += 1
                print(count_N1,count_N2)
            
                                                                                           
        if N == 2:
            #print(2)
            F = N + S
                
            #Geometry
                
            for distance in distance_list:                    
                k0=2*np.pi
                k_vec = [0,k0,0]
                k = np.linalg.norm(k_vec)
                k_unity = np.divide(k_vec,k)
                n_hat = np.divide(n_hat,np.linalg.norm(n_hat))
                k_vec = k0*k_unity
                d = distance/k0
                p = [0,0,1]
                P = [p]*N
                R = position_vector.position_vector(d,N)
                
                for Omega in Omega_list:
                       
                    path_data = "../results/data/QT/Matrix_hdf12_Correlation_Intensity_Conditional_N=%d_kd=%.2f_Omega=%.2f_Delta=%.2f/"%(N,k0*d,Omega,Delta)
                        
                    filename = path_data+filename_generator.filename_gen_sensor("I",N,S,k,d,Omega,Delta)
                    f = h5py.File(filename,'r')
                    g2_Matrix = f['/1/g2']

                    R = np.zeros((omega_points,omega_points))
                    conditional_Matrix = np.zeros((omega_points,omega_points))
                    count = 0
                    count_list = []
    
                    for i in range(0,omega_points):
                        for j in range(0,omega_points):
                            R[i,j] = ((g2_Matrix[i,j])**(2))/(g2_Matrix[i,i]*g2_Matrix[j,j])

                    g2_Condition = np.array(list(map(lambda x: x>10,g2_Matrix)))

                    R_Condition = np.array(list(map(lambda y: y>1,R)))
                    Conditional_Matrix = np.zeros((omega_points,omega_points))
                    M = np.zeros((omega_points,omega_points))

                    for i in range(0,omega_points):
                        for j in range(0,omega_points):
                            Conditional_Matrix[i,j] = g2_Condition[i,j] and R_Condition[i,j]
                            if Conditional_Matrix[i,j] == True:
                                M[i,j] = g2_Matrix[i,j]
                            else: 
                                M[i,j] = float('nan')
                
                    vmin,vmax = vmin_g,vmax_g
                    x = omega_list
                    y = omega_list
                    title = 'Conditional Map - CSI'
                    x_title = '$\omega_{2}$'
                    y_title = '$\omega_{1}$'
                    file_name_map = 'Conditional_map_CSI_N=%d_kd=%.2f_Omega=%.2f_Delta=%.2f.pdf'%(N,k0*d,Omega,Delta)
                    file_name = path_plot + file_name_map
                    Two_Dim_Maps.TwoDim_MapLogScale(M,N,k0,d,x,y,x_title,y_title,file_name,title,vmin,vmax)            
                    count_N2 += 1
                    print(count_N1,count_N2)
                    
                        
                        
                        

                    
                                                                                                                                
                













