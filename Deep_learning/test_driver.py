
# ================== ABOUT ====================
# Author: Jean Bonnemain
# Modifications: Matthias Zeller

# Information: run the script with '-h' or help to print informations !

# ================= IMPORTS ===================

from __future__ import absolute_import, division, print_function, unicode_literals

# TensorFlow and tf.keras
import tensorflow as tf
from tensorflow import keras

# Helper libraries
import numpy as np
import matplotlib.pyplot as plt

import scipy
import scipy.io as sio

import sys
sys.path.insert(1, sys.path[0] + '/../Simulation_script/')

import utils_openmodelica as uo
import utils_deeplearning as udl
import datetime
import os

# ============== PARAMETERS ===============

# PERCENTAGE OF COEFFICIENTS TO KEEP
perccoef = 0.05

# ============ IMPLEMENTATION =============

# One can run the script with arguments to run specific parts
run_train, run_test, path = udl.manage_args(sys.argv)

print('Tensorflow version:', tf.__version__)

# If we only want to train the DNN, we also want to store the test data
train_only = True if run_train and not run_test else False
# If we only want to test the DNN, we have to load the the test data
test_only  = True if run_test and not run_train else False

# Training mode
if run_train:
    print("\n========== TRAINING")
    # Launch DNN training
    res = udl.train_dnn(perccoef, files_path=path, save_test_data=train_only)

    # If we train & test, load data from memory rather than from files
    if not train_only:
        normdata = res[0]
        (Xtest, Ytest) = res[1]
        model = res[2]
        del res

# Testing mode
if run_test:
    # If run-only mode, load data from files
    if test_only:
        try:
            Xtest = np.load('Xtest.npy')
            Ytest = np.load('Ytest.npy')
        except FileNotFoundError:
            print("\nERROR ! Xtest or Ytest not found, you probably did not run the training mode before")
            exit()

        normdata = {
            'coefmins':None, 'coefmaxs':None,
            'parammins':None, 'parammaxs':None
        }
        for name in normdata:
            normdata[name] = np.load(name + '.npy')

        model = keras.models.load_model("DNN_0D_Model.h5")
        print("Data were loaded from files...")
        for name,val in normdata.items():
            print(name, '=', val)
        print("Xtest.shape =", Xtest.shape)
        print("Ytest.shape =", Ytest.shape)


    print("\n========== TESTING")

    # Launc DNN testing
    udl.test_dnn(model, Xtest, Ytest, normdata)


print("\nExiting to keep control...")
exit()

filepath="/Users/jean.bonnemain/Documents/Code/0d_model/Modelica_Code/0D_Original/"
today = datetime.datetime.now()
outputfolder=os.getcwd() + '/' + today.strftime("%Y") + '_' + today.strftime("%m") + '_' + today.strftime("%d")

uo.prepareOutputFolder(outputfolder)

testsamples = Ypred.shape[0]

# run 0D simulations with OpenModelica to compare results of the 0D with exact
# values and predicted values

for indexsample in range(0,testsamples):
    param1 = uo.Parameter('Param_LeftVentricle_Emax0')
    param1.setValue(Ypred[indexsample,0])
    param2 = uo.Parameter('Param_LeftVentricle_EmaxRef0')
    param2.setValue(Ypred[indexsample,1])
    param3 = uo.Parameter('Param_LeftVentricle_AGain_Emax')
    param3.setValue(Ypred[indexsample,2])
    param4 = uo.Parameter('Param_LeftVentricle_kE')
    param4.setValue(Ypred[indexsample,3])

    listParameters = [param1, param2, param3, param4]

    suffix = str(indexsample) + "predicted"
    uo.launchSimulation(filepath, listParameters, suffix, outputfolder)
    
    param1 = uo.Parameter('Param_LeftVentricle_Emax0')
    param1.setValue(Ytest[indexsample,0])
    param2 = uo.Parameter('Param_LeftVentricle_EmaxRef0')
    param2.setValue(Ytest[indexsample,1])
    param3 = uo.Parameter('Param_LeftVentricle_AGain_Emax')
    param3.setValue(Ytest[indexsample,2])
    param4 = uo.Parameter('Param_LeftVentricle_kE')
    param4.setValue(Ytest[indexsample,3])

    listParameters = [param1, param2, param3, param4]

    suffix = str(indexsample) + "exact"
    uo.launchSimulation(filepath, listParameters, suffix, outputfolder)