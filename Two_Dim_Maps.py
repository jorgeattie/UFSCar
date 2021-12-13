# -*- coding: utf-8 -*-
"""
Created on Fri Mar 19 23:56:01 2021

@author: Jorge
"""

from qutip import *
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plot
import scipy.linalg as la
import matplotlib.pyplot as plt

from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
from matplotlib import ticker, cm
import matplotlib.colors as colors
import matplotlib.cbook as cbook






def TwoDim_Map(M,N,k0,d,x,y,x_title,y_title):
    #nrm = mpl.colors.Normalize(0, 2)
    fig, axes = plot.subplots(1,1,figsize=(6, 3))
    axes.set_title("$N = %d \quad kd = %.2f$" %(N,k0*d));
    axes.set_xlabel(x_title)
    axes.set_ylabel(y_title)
    plt = axes.contourf(x,y,M,100, cmap=cm.gist_rainbow)
    plt.rcParams.update({'font.size':14})
    cb1 = fig.colorbar(plt, ax=axes)


def TwoDim_MapLogScale(M,N,k0,d,x,y,x_title,y_title,file_name,title):
    fig, axes = plot.subplots(1,1,figsize=(6,3))
    axes.labelsize:xx-large
    #axes.set_title("$N = %d \quad kd = %.2f$" %(N,k0*d));
    axes.tick_params('x',labelsize=16)
    axes.tick_params('y',labelsize=16)
    axes.set_xlabel(x_title,fontsize=18)
    axes.set_ylabel(y_title,fontsize=18)       
    pcm = axes.pcolor(x,y,M,shading='auto',norm=colors.LogNorm(vmin=10**(-21),vmax=10**(-10)),cmap=cm.gist_rainbow)
    #pcm = axes.pcolor(x,y,M,shading='auto',cmap=cm.rainbow)
    #fig.colorbar(pcm,ax=axes,extend='max',label = title)
    fig.colorbar(pcm,ax=axes,extend='max',label=title)
    plt.rcParams.update({'font.size':14})
    plt.savefig(file_name,bbox_inches = 'tight')
    #plt.show()
    


#def Two_Dim_Maps(x,y,M,N,k0,d,x_title,y_title):        
#    nrm = mpl.colors.Normalize(0, 2)
#    fig, axes = plot.subplots(1,1,figsize=(10, 4))
#    #plt = axes.contourf(x,y,M,100, cmap=cm.RdBu, norm=nrm)
#    plt = axes.contourf(x,y,M,100, cmap=ListedColormap(["darkorange", "gold", "lawngreen", "lightseagreen"]) , norm=nrm)
#    axes.set_title("$N = %d \quad kd = %.2f$" %(N,k0*d));
#    axes.set_xlabel(x_title)
#    axes.set_ylabel(y_title)
#    cb1 = fig.colorbar(plt, ax=axes)
