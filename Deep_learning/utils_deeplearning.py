#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 21 15:11:02 2020

@author: jean.bonnemain
"""

import numpy as np

# Normalizes output parameters. If we want to use max and min parameters from 
# the trained model, just give them as input.
def normalizeoutputmatDL(mat,newmins,newmaxs, parammins=None, parammaxs=None):
    shapes = mat.shape
    
    if type(parammins) == type(None) and type(parammaxs) == type(None):
        parammins = np.full(newmins.shape,0.0)
        parammaxs = np.full(newmins.shape,0.0)
        for i in range(0,shapes[1]):
            param = mat[:,i]
            parammin = param.min()
            parammax = param.max()
            parammins[i] = parammin
            parammaxs[i] = parammax
            mat[:,i] = newmins[i] + (newmaxs[i] - newmins[i]) * (param - parammin) / (parammax - parammin)
            print("Compute min and max for output parameters")

    else:
        for i in range(0,shapes[1]):
            param = mat[:,i]
            parammin = parammins[i]
            parammax = parammaxs[i]
            mat[:,i] = newmins[i] + (newmaxs[i] - newmins[i]) * (param - parammin) / (parammax - parammin)
            print("load min and max for output parameters")

    return mat, parammins, parammaxs

# Normalizes output parameters. If we want to use max and min parameters from 
# the trained model, just give them as input.
def normalizeinputmatDL(mat,newmins,newmaxs,coefmins=None,coefmaxs=None):
    shapes = mat.shape
    
    if type(coefmins) == type(None) and type(coefmaxs) == type(None):
        coefmins = np.full(newmins.shape,0.0)
        coefmaxs = np.full(newmaxs.shape,0.0)
        for coefind in range(0,shapes[1]):
            for freq in range(0,shapes[2]):
                coef = mat[:,coefind,freq]
                coefmin = coef.min()
                coefmax = coef.max()
                coefmins[coefind,freq] = coefmin
                coefmaxs[coefind,freq] = coefmax
                if (abs(coefmax - coefmin) > 1e-15):#In case of normalization by column (bk removed in MATLAB)
                    mat[:,coefind,freq] = newmins[coefind,freq] + (newmaxs[coefind,freq] - newmins[coefind,freq]) * (coef - coefmin) / (coefmax - coefmin)
                    print("Compute min and max for input parameters")

    else:
        for coefind in range(0,shapes[1]):
            for freq in range(0,shapes[2]):
                coef = mat[:,coefind,freq]
                coefmin = coefmins[coefind,freq]
                coefmax = coefmaxs[coefind,freq]
                if (abs(coefmax - coefmin) > 1e-15):#In case of normalization by column (bk removed in MATLAB)
                    mat[:,coefind,freq] = newmins[coefind,freq] + (newmaxs[coefind,freq] - newmins[coefind,freq]) * (coef - coefmin) / (coefmax - coefmin)
                    print("Load min and max for input parameters")
    return mat, coefmins, coefmaxs