#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 20 15:59:13 2025

@author: baptiste
"""
import numpy as np 

#return j, a, b

def getCoeff(num_left, num_right, num_down, num_up, num_cent, type_cent, cl_cent):
    if type_cent == 0:
        raise ValueError("Invalid value for type_cent: 0 is not in the domain of calcul.")
    
    if type_cent == 2:
        a = np.array([1]).reshape(-1,1)
        j = np.array([num_cent]).reshape(-1,1)
        b = cl_cent
    
    if type_cent == 1:
        j = np.array([num_left, num_right, num_down, num_up, num_cent]).reshape(-1,1)  # Transpose the array
        a = np.array([1, 1, 1, 1, -4]).reshape(-1,1)
        b = 0
    
    return j , a , b