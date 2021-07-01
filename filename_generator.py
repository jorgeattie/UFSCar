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

def filename_gen(TYPE,N,k0,d):
    return "%s_N=%d_kd=%.2f.h5" %(TYPE,N,k0*d)