#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 20 15:59:13 2025

@author: baptiste
"""
import numpy as np 
import scipy as sc

#return j, a, b

def getCoeff(num_left, num_right, num_down, num_up, num_cent, type_cent, cl_cent) :
    if(type_cent == 0) :
        raise ValueError("Invalid value for type_cent: 0 is not in the domain of calcul.")
        return 0, 0, 0
    if(type_cent == 2) :
        a = np.array(1)
        j = np.array(num_cent)
        return j, a, cl_cent
    if(type_cent == 1) :
        j = np.array([num_left, num_right, num_down, num_up, num_cent])
        a = np.array([1, 1, 1, 1, -4])
        return j, a, cl_cent
 