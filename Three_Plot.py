# -*- coding: utf-8 -*-
"""
Created on Mon Mar  8 14:28:31 2021

@author: Jorge
"""

from qutip import *
import numpy as np
import matplotlib.pyplot as plt
import scipy.linalg as la

def Three_Plot(x,y1,y2,y3,color1,color2,color3,title1,title2,title3,x_title,y_title,scale,N,k0,d,Omega,rnge,name):
    fig,axes = plt.subplots(1,1,figsize=(12,6))
    axes.plot(x,y1,color1,label = title1)
    axes.plot(x,y2,color2,label = title2)
    axes.plot(x,y3,color3,label = title3)
    plt.ylim(rnge)
    axes.legend(loc=0)
    axes.set_xlabel(x_title)
    axes.set_ylabel(y_title)
    plt.yscale(scale)
    plt.title(name)
    plt.show()