from qutip import *
import numpy as np
import matplotlib.pyplot as plt
import scipy.linalg as la

def Plot(x,y,color,title,x_title,y_title,scale,N,k0,name,rnge_x,rnge_y,file_name):
    fig,axes = plt.subplots(1,1,figsize=(6,3))
    #axes.plot(x,y,color,label = title)
    axes.plot(x,y,color)
    plt.xlim(rnge_x)
    plt.ylim(rnge_y)
    axes.legend(loc=0)
    axes.set_xlabel(x_title)
    axes.set_ylabel(y_title)
    plt.yscale(scale)
    plt.title(name)
    plt.rcParams.update({'font.size':12})
    #plt.show()
    plt.savefig(file_name,bbox_inches = 'tight')

def Plot_wHorizontalLine(x,y,color,x_title,y_title,scale,N,k0,d,name,rnge_x,rnge_y,y_ref,file_name):
    fig,axes = plt.subplots(1,1,figsize=(6,3))
    plt.xticks(fontsize=14)
    axes.plot(x,y,color)
    #axes.axis('equal')
    r_x = [np.min(rnge_x),np.max(rnge_x)]
    if N == 3:
        r_y = [10**(-6),np.max(rnge_y)*10]
    
    if N == 2:
        r_y = [10**(-6),10**(4)]
        
    plt.xlim(r_x)
    plt.ylim(r_y)
    #labels_x = rnge_x 
    #labels_y = [1,2,3]
    #plt.xticks(np.arange(np.min(rnge_x),np.max(rnge_x),500))
    plt.yticks(np.array([10**(-3),10**(0),10**(3)]))
    axes.tick_params('x',labelsize=16)
    axes.tick_params('y',labelsize=16)
    axes.set_xlabel(x_title,fontsize=18)
    axes.set_ylabel(y_title,fontsize=18)
    plt.yscale(scale)
    plt.axhline(y_ref,xmin=np.min(x),xmax=np.max(x),color='black',linestyle='--')
    #plt.title(name)
    plt.rcParams.update({'font.size':14})
    plt.savefig(file_name,bbox_inches = 'tight')
    #plt.show()
    
def Plot_GCurves(x,y,color,x_title,y_title,scale,name):
    fig,axes = plt.subplots(1,1,figsize=(8,4))
    axes.plot(x,y,color)
    #axes.legend(loc=0)
    axes.set_xlabel(x_title)
    axes.set_ylabel(y_title)
    plt.yscale(scale)
    plt.title(name)
    plt.rcParams.update({'font.size':22})
    plt.show()
    