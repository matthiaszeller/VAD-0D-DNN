
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

import utils_deeplearning as udl
import utils_openmodelica as uo
from setup import *
import datetime
import os

# ============== PARAMETERS ===============

# PERCENTAGE OF COEFFICIENTS TO KEEP
perccoef = 0.05

# ============ IMPLEMENTATION =============

# Manage arguments
args = udl.manage_args(sys.argv[1:])

print('Tensorflow version:', tf.__version__)

# Training mode
if args.mode == 'train':
    print("\n========== TRAINING")

    # Launch DNN training
    try:
        res = udl.train_dnn(perccoef, files_path=args.path, save_test_data=True,
                            verbose=args.verbose, selected_aks=args.selectedaks)
    except FileNotFoundError:
        print('\n<ERROR> X.mat or Y.mat was not found, you should use the --path argument.')


# Testing mode
elif args.mode == 'test':
    # Load data from files
    print("\n========== LOADING")
    try:
        # dnn_folder defined in Simulation_script/setup.py
        Xtest = np.load(dnn_folder + '/Xtest_norm.npy')
        Ytest = np.load(dnn_folder + '/Ytest_norm.npy')
    except FileNotFoundError:
        print("\nERROR ! Xtest_norm and/or Ytest_norm not found, you probably \
        did not run the training mode before or 'dnn_folder' in Simulation_script/setup.py does not target the right path")
        exit()

    normdata = {
        'coefmins':None, 'coefmaxs':None,
        'parammins':None, 'parammaxs':None
    }
    for name in normdata:
        normdata[name] = np.load(dnn_folder + '/' + name + '.npy')

    model = keras.models.load_model(dnn_folder + "/DNN_0D_Model.h5")

    print("Data were loaded from files...")

    for name,val in normdata.items():
        print(name, '=', val)
    print("Xtest.shape =", Xtest.shape)
    print("Ytest.shape =", Ytest.shape)


    print("\n========== TESTING")

    # Parameters of the simulation
    # WARNING: the order of the parameters defined below must be consistent with
    # the order defined in the simulation script (launch_openmodelica.py)
    param1 = uo.Parameter('Param_LeftVentricle_Emax0')
    param2 = uo.Parameter('Param_LeftVentricle_EmaxRef0')
    param3 = uo.Parameter('Param_LeftVentricle_AGain_Emax')
    param4 = uo.Parameter('Param_LeftVentricle_kE')
    param_lst = [param1, param2, param3, param4]

    # Setup the path for the model outputs
    today = datetime.datetime.now()
    # output_folder_DNN_test defined in Simulation_script/setup.py
    output_folder_DNN_test += '_' + today.strftime("%Y") + '_' + today.strftime("%m") + '_' + today.strftime("%d")

    # Launch DNN testing
    udl.test_dnn(model, Xtest, Ytest, normdata, param_lst,
                 output_folder_DNN_test, file_path, dnn_folder)

else:
    print('Error')


