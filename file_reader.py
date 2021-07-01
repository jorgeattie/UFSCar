# -*- coding: utf-8 -*-
"""
Created on Fri Apr 30 16:52:16 2021

@author: Jorge
"""

from qutip import *
import numpy as np
import matplotlib.pyplot as plt
import scipy.linalg as la

import h5py 
import os 

import filename_generator

filename_QT = path_data_QT+filename_generator.filename_gen('QT', N, k0, d)
f = h5py.File(filename_QT,'r')
M = f['/1/M']
t_list = f['/1/t']

M_QT = Qobj(M[t],dims=[[2 for i in range(N)],[2 for i in range(N)]])