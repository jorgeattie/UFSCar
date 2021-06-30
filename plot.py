# -*- coding: utf-8 -*-
"""
Created on Fri Mar  5 16:14:13 2021

@author: Jorge
"""

from qutip import *
import numpy as np
import matplotlib.pyplot as plt
import scipy.linalg as la

def Plot(x,y,color,title,x_title,y_title,scale,N,k0,d,name,rnge):
    fig,axes = plt.subplots(1,1,figsize=(12,6))
    axes.plot(x,y,color,label = title)
    plt.ylim(rnge)
    axes.legend(loc=0)
    axes.set_xlabel(x_title)
    axes.set_ylabel(y_title)
    plt.yscale(scale)
    plt.title(name)
    plt.show()
    