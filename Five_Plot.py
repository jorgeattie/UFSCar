# -*- coding: utf-8 -*-
"""
Created on Mon Mar  8 16:38:55 2021

@author: Jorge
"""

from qutip import *
import numpy as np
import matplotlib.pyplot as plt
import scipy.linalg as la

def Five_Plot(x,y1,y2,y3,y4,y5,color1,color2,color3,color4,color5,title1,title2,title3,title4,title5,x_title,y_title,scale,N,k0,d,Omega,rnge):
    fig,axes = plt.subplots(1,1,figsize=(12,6))
    axes.plot(x,y1,color1,label = title1)
    axes.plot(x,y2,color2,label = title2)
    axes.plot(x,y3,color3,label = title3)
    axes.plot(x,y4,color4,label = title4)
    axes.plot(x,y5,color5,label = title5)
    plt.ylim(rnge)
    axes.legend(loc=0)
    axes.set_xlabel(x_title)
    axes.set_ylabel(y_title)
    plt.yscale(scale)
    plt.title("$N = %d \quad kd = %.2f \quad \Omega = %.1f$" %(N,k0*d,Omega))
    plt.show()