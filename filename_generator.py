# -*- coding: utf-8 -*-
"""
Created on Fri Apr 30 16:49:54 2021

@author: Jorge
"""

from qutip import *
import numpy as np
import matplotlib.pyplot as plt
import scipy.linalg as la

import h5py 
import os 

def filename_gen_sensor(TYPE,N,S,k,d,Omega,Delta):
    return "%s_N=%d_S=%d_kd=%.2f_Omega=%.2f_Delta=%.2f.h5" %(TYPE,N,S,k*d,Omega,Delta)

def filename_gen_pop_ratio(TYPE,N,k,d):
    return "%s_N=%d_kd=%.2f"%(TYPE,N,k*d)