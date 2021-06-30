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

from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap




def TwoDim_Map(M,N,k0,d,x,y,x_title,y_title):
    #nrm = mpl.colors.Normalize(0, 2)
    fig, axes = plot.subplots(1,1,figsize=(10, 4))
    axes.set_title("$N = %d \quad kd = %.2f$" %(N,k0*d));
    axes.set_xlabel(x_title)
    axes.set_ylabel(y_title)
    plt = axes.contourf(x,y,M,100, cmap=cm.gist_rainbow)
    cb1 = fig.colorbar(plt, ax=axes)

#def Two_Dim_Maps(x,y,M,N,k0,d,x_title,y_title):        
#    nrm = mpl.colors.Normalize(0, 2)
#    fig, axes = plot.subplots(1,1,figsize=(10, 4))
#    #plt = axes.contourf(x,y,M,100, cmap=cm.RdBu, norm=nrm)
#    plt = axes.contourf(x,y,M,100, cmap=ListedColormap(["darkorange", "gold", "lawngreen", "lightseagreen"]) , norm=nrm)
#    axes.set_title("$N = %d \quad kd = %.2f$" %(N,k0*d));
#    axes.set_xlabel(x_title)
#    axes.set_ylabel(y_title)
#    cb1 = fig.colorbar(plt, ax=axes)
