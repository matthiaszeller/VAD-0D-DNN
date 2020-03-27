#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 21 15:11:02 2020

@author: jean.bonnemain
"""

import numpy as np
import scipy.io as sio
import os.path

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


def train_dnn(perccoef, files_path=''):
    # ======== DATA LOADING
    X = sio.loadmat(files_path+'X.mat')['X']
    Y = sio.loadmat(files_path+'Y.mat')['Y']
    print("Loaded X.mat ({}) and Y.mat ({}) files"
          .format('x'.join([str(n) for n in X.shape]),
                  'x'.join([str(n) for n in Y.shape])))

    ncoefficients = X.shape[1]
    nfrequencies = X.shape[2]
    noutparams = Y.shape[1]
    nsamples = X.shape[0]

    # ======== DATA REDUCTION
    # Select a subset of frequencies for X according to argument 'perccoef'
    if (perccoef < 0.99999):
        if nfrequencies % 2 == 0:
            aks = (nfrequencies + 2) / 2
            bks = aks - 2
        else:
            aks = (nfrequencies + 1) / 2
            bks = aks - 1
        aks = int(aks)
        selectedaks = int(np.floor(aks * perccoef))
        indicestoselect = list(range(0, selectedaks)) + list(range(aks, aks + selectedaks - 1))
        X = np.take(X, indicestoselect, axis=2)
        nfrequencies = X.shape[2]
        print("Reduced amount of input, X is now " + str(X.shape))

    # ======== NORMALIZATION
    newmins = np.full([ncoefficients, nfrequencies], 0.0)
    newmaxs = np.full([ncoefficients, nfrequencies], 1.0)

    X, coefmins, coefmaxs = udl.normalizeinputmatDL(X, newmins, newmaxs)
    np.save("coefmins", coefmins)
    np.save("coefmaxs", coefmaxs)

    newmins = np.full([noutparams], 0.0)
    newmaxs = np.full([noutparams], 1.0)
    Y, parammins, parammaxs = udl.normalizeoutputmatDL(Y, newmins, newmaxs)
    np.save("parammins", parammins)
    np.save("parammaxs", parammaxs)



def manage_args(args):
    N = len(args)
    # Running mode
    run_train, run_test = True, True
    # Location of X.mat and Y.mat
    path_files = ''

    # Get rid of the python script in the args
    args = args[1:]

    # Specification of the files path
    if '--path' in args:
        i = args.index('--path')
        folder = args[i+1]
        if folder[-1] != '/': folder += '/'
        if not os.path.isdir(folder):
            print("ERROR ! Invalid directory '{}'".format(folder))
            exit()
        path_files = folder

    # Help
    if '-h' in args or 'help' in args:
        print("The following commands can be provided (all are optional, if none is provided -> test & train by default):")
        #print("* 'help' or '-h'")
        print("* 'train'\t\tTrain the DNN with X.mat, Y.mat.\n\t\t\tCreates files 1-6.")
        print("* 'test'\t\tTest the DNN by loading files X, X, X")
        print("* '--path <folder>'\tSpecify the directory of X.mat and Y.mat")
        exit()
    # Train only
    if 'train' in args:
        run_test = False
        print("Running in train mode...")
    # Test only
    elif 'test' in args:
        run_train = False
        print("Running in test mode...")
    # Error
    #else:
    #    print("ERROR ! Invalid argument(s).\nProgram abortion.")
    #    exit()

    if run_test and run_train:
        print("Running train and test mode...")

    return run_train, run_test, path_files