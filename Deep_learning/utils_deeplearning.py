#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 21 15:11:02 2020

@author: jean.bonnemain
"""

import numpy as np
import scipy.io as sio
import os.path
from tensorflow import keras
import matplotlib.pyplot as plt
import utils_openmodelica as uo

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


# for simplicity this is specific to a three dimensional tensor in which the
# the second component differentiates between aks and bks
def normalizeinputmat(mat, newmins, newmaxs):
    #    # compute min and max values of aks
    #    mataks = mat[:,0,:]
    #    aksmin = mataks.min()
    #    aksmax = mataks.max()
    #
    #    # compute min and max values of bks
    #    matbks = mat[:,1,:]
    #    bksmin = matbks.min()
    #    bksmax = matbks.max()
    #
    #    # normalize such that aks and bks are in range (newmin,newmax)
    #    mat[:,0,:] = newmin + (newmax - newmin) * (mataks - aksmin) / (aksmax - aksmin)
    #    mat[:,1,:] = newmin + (newmax - newmin) * (matbks - bksmin) / (bksmax - bksmin)
    shapes = mat.shape

    coefmins = np.full(newmins.shape, 0.0)
    coefmaxs = np.full(newmaxs.shape, 0.0)

    for coefind in range(0, shapes[1]):
        for freq in range(0, shapes[2]):
            coef = mat[:, coefind, freq]
            coefmin = coef.min()
            coefmax = coef.max()
            coefmins[coefind, freq] = coefmin
            coefmaxs[coefind, freq] = coefmax
            if (abs(coefmax - coefmin) > 1e-15):  # In case of normalization by column (bk removed in MATLAB)
                mat[:, coefind, freq] = newmins[coefind, freq] + (newmaxs[coefind, freq] - newmins[coefind, freq]) * (
                            coef - coefmin) / (coefmax - coefmin)

    return mat, coefmins, coefmaxs


# different normalization: we normalize wrt to the min and max values of each
# parameter at different samples
def normalizeoutputmat(mat, newmins, newmaxs):
    shapes = mat.shape
    parammins = np.full(newmins.shape, 0.0)
    parammaxs = np.full(newmins.shape, 0.0)

    for i in range(0, shapes[1]):
        param = mat[:, i]
        parammin = param.min()
        parammax = param.max()
        parammins[i] = parammin
        parammaxs[i] = parammax
        mat[:, i] = newmins[i] + (newmaxs[i] - newmins[i]) * (param - parammin) / (parammax - parammin)

    return mat, parammins, parammaxs


def build_keras_model(input_shape, noutparams):
    model = keras.Sequential([
        keras.layers.Flatten(input_shape=input_shape),
        keras.layers.Dense(32, activation='relu'),
        keras.layers.Dense(32, activation='relu'),
        keras.layers.Dense(32, activation='relu'),
        keras.layers.Dense(32, activation='relu'),
        keras.layers.Dense(32, activation='relu'),
        keras.layers.Dense(noutparams, activation='sigmoid')
    ])

    model.compile(optimizer='adam',
                  loss='mse',
                  metrics=['mae'])

    return model

def create_and_save_performance_fig(Ytest, Ypred, normdata, folder):
    fig, axs = plt.subplots(2, 2)
    axs[0, 0].scatter(Ytest[:, 0], Ypred[:, 0])
    axs[0, 0].set_title('Left Ventricle Emax0')
    axs[0, 0].set_xlabel('real parameter')
    axs[0, 0].set_ylabel('predicted parameter')
    axs[0, 1].scatter(Ytest[:, 1], Ypred[:, 1])
    axs[0, 1].set_title('Left Ventricle EmaxRef0')
    axs[0, 1].set_xlabel('real parameter')
    axs[0, 1].set_ylabel('predicted parameter')
    axs[1, 0].scatter(Ytest[:, 2], Ypred[:, 2])
    axs[1, 0].set_title('Left Ventricle AGain_Emax')
    axs[1, 0].set_xlabel('real parameter')
    axs[1, 0].set_ylabel('predicted parameter')
    axs[1, 1].scatter(Ytest[:, 3], Ypred[:, 3])
    axs[1, 1].set_title('Left Ventricle kE')
    axs[1, 1].set_xlabel('real parameter')
    axs[1, 1].set_ylabel('predicted parameter')
    axs[1, 1].set_xlim(normdata['parammins'][3], normdata['parammaxs'][3])
    axs[1, 1].set_ylim(normdata['parammins'][3], normdata['parammaxs'][3])
    fig.set_figheight(10)
    fig.set_figwidth(15)
    plt.savefig(folder + '/DNN_Performance.eps')

def test_dnn(model, Xtest, Ytest, normdata, param_lst,
             output_dnn_test, modelica_file_path, dnn_folder):
    # Xtest, Ytest are normalizedf

    # ======== PREDICTION
    Ypred = model.predict(Xtest)

    # ======== ORIGINAL DATA RECOVERY
    # "Un-normalize" values
    Xtest, _, _ = normalizeinputmat(Xtest, normdata['coefmins'], normdata['coefmaxs'])
    Ypred, _, _ = normalizeoutputmat(Ypred, normdata['parammins'], normdata['parammaxs'])
    Ytest, _, _ = normalizeoutputmat(Ytest, normdata['parammins'], normdata['parammaxs'])
    print(Ypred[0, :])
    print(Xtest[0, 0, :])  # return frequency for systemic arteries

    # scipy.io.savemat('Ypred.mat',mdict={'Ypred': Ypred})
    # scipy.io.savemat('Ypred.mat',mdict={'Xtest': Xtest})

    # ======== GENERATE AND SAVE FIGURE OF DNN PERFORMANCE
    create_and_save_performance_fig(Ytest, Ypred, normdata, dnn_folder)
    print("exit to keep control")
    exit()
    # ======== SIMULATE WITH MODELICA
    uo.runTestSimulation(Ytest=Ytest, Ytest_pred=Ypred, param_lst=param_lst,
                         output_dnn_test=output_dnn_test, modelica_file_path=modelica_file_path)

def train_dnn(perccoef, files_path='', save_test_data=True):
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

    # ======== (INPUT) DATA REDUCTION
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

    # cut the input matrix (e.g. to keep only a subset of the frequencies)
    # freqmax = 100
    # X = X[:,:,0:freqmax]

    # ======== NORMALIZATION
    newmins = np.full([ncoefficients, nfrequencies], 0.0)
    newmaxs = np.full([ncoefficients, nfrequencies], 1.0)

    X, coefmins, coefmaxs = normalizeinputmatDL(X, newmins, newmaxs)

    newmins = np.full([noutparams], 0.0)
    newmaxs = np.full([noutparams], 1.0)
    Y, parammins, parammaxs = normalizeoutputmatDL(Y, newmins, newmaxs)

    # Save normalization values
    normdata = {
        'coefmins':coefmins,
        'coefmaxs':coefmaxs,
        'parammins':parammins,
        'parammaxs':parammaxs
    }
    for name, data_elem in normdata.items():
        np.save(name, data_elem)

    # ======== DATA SPLITTING: TRAINING, VALIDATION, TEST SETS
    perctest = 0.05
    perctrain = 1 - perctest
    percvalidation = 0.2 # percentage applied to perctrain and not to the 1.00
    samplestrain = int(np.floor(perctrain * nsamples))
    samplestest = nsamples - samplestrain

    Xtrain = X[0:samplestrain,:,:]
    Ytrain = Y[0:samplestrain,:]
    Xtest = X[samplestrain+1:,:,:]
    Ytest = Y[samplestrain+1:,:]

    data = {'Xtrain':Xtrain, 'Xtest':Xtest, 'Ytrain':Ytrain, 'Ytest':Ytest}
    print("Training and test set generated: shapes are", {key:val.shape for key,val in data.items()})

    # ======== BUILD KERAS MODEL
    model = build_keras_model(input_shape=(ncoefficients, nfrequencies), noutparams=noutparams)
    print(model.summary())

    # ======== TRAIN DNN
    print("Fitting the model...")
    history = model.fit(Xtrain, Ytrain, epochs=1000, validation_split=percvalidation, verbose=0)
    print("Model is trained")

    # ======== SAVE DATA
    # Create and save figure of loss and validation loss
    plt.semilogy(history.history['loss'], label='loss')
    plt.semilogy(history.history['val_loss'], label='val_loss')
    plt.legend()
    plt.savefig('Losses.eps')

    # Save the model in the .h5 format
    model.save('DNN_0D_Model.h5')

    # If train-only mode, save X and Y into files and return nothing (training and testing sets)
    if save_test_data:
        np.save('Xtest_norm', Xtest)
        np.save('Ytest_norm', Ytest)

    #  Otherwise, do not save data but return them from memory
    return (normdata, (Xtest, Ytest), model)


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
        print("\nThe following commands can be provided (all are optional, if none is provided -> test & train by default):")
        #print("* 'help' or '-h'")
        print("* 'train'\t\tTrain the DNN with X.mat, Y.mat.\n\t\t\tCreates files 1-8.")
        print("* 'test'\t\tTest the DNN by loading files 1-4 and 6-8")
        print("* '--path <folder>'\tSpecify the directory of X.mat and Y.mat")
        print("\nFiles description:")
        print("1. coefmins\tminimums of the original data (input) before normalization")
        print("2. coefmaxs\tmaximums of the original data (input) before normalization")
        print("3. parammins\tminimums of the original data (output) before normalization")
        print("4. parammaxs\tmaximums of the original data (output) before normalization")
        print("5. Losses.eps\tfigure showing the MSE of the training and validation sets")
        print("6. DNN_0D_Model.h5\tDNN model")
        print("7. Xtest.npy\tInput test data")
        print("8. Ytest.npy\tOutput test data")
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