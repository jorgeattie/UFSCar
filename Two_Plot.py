# -*- coding: utf-8 -*-
"""
Created on Mon Mar  8 14:20:10 2021

@author: Jorge
"""

from qutip import *
import numpy as np
import matplotlib.pyplot as plt
import scipy.linalg as la

def Two_Plot(x,y1,y2,color1,color2,title1,title2,x_title,y_title,scale,N,k0,d,name,rnge):
    fig,axes = plt.subplots(1,1,figsize=(12,6))
    axes.plot(x,y1,color1,label = title1)
    axes.plot(x,y2,color2,label = title2)
    plt.ylim(rnge)
    axes.legend(loc=0)
    axes.set_xlabel(x_title)
    axes.set_ylabel(y_title)
    plt.yscale(scale)
    plt.title(name)
    plt.show()