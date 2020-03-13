from __future__ import absolute_import, division, print_function, unicode_literals

# TensorFlow and tf.keras
import tensorflow as tf
from tensorflow import keras

# Helper libraries
import numpy as np
import matplotlib.pyplot as plt

import scipy
import scipy.io as sio

import matplotlib.pyplot as plt

import sys 
sys.path.insert(1, '../Simulation_script/')

import utils_openmodelica as uo
import utils_deeplearning as udl

import datetime
import os

print('Version Tensorflow:')
print(tf.__version__)

modelfile="/Users/jean.bonnemain/switchdrive/Institution/Results/2020_01_21_Test_Procedure_load_DNN/Train_set/DNN_0D_Model.h5"
#filepath="/Users/jean.bonnemain/switchdrive/Institution/Results/2019_12_24/"
#filepath="/Users/jean.bonnemain/switchdrive/Institution/Results/2020_01_17_Heart_failure_parameters/2020_01_14_4_params/"
#filepath="/Users/jean.bonnemain/switchdrive/Institution/Results/2020_01_17_Heart_failure_parameters/2020_01_14_4_params/100_run/"

# Load the DNN model from .h file
dnn_model = tf.keras.models.load_model(modelfile)

# Show the model architecture
dnn_model.summary()

# Load data
X = sio.loadmat('X.mat')['X']
Y = sio.loadmat('Y.mat')['Y']

# Load min and max coefficients for normalization (same as the trained model)
parammins = np.load("parammins.npy")
parammaxs = np.load("parammaxs.npy")
coefmins = np.load("coefmins.npy")
coefmaxs = np.load("coefmaxs.npy")


# Select a subset of frequencies
perccoef = 0.05
ncoefficients = X.shape[1] # Number of parameters (arterial and pulmonary curves)
nfrequencies  = X.shape[2] # Fourier coefficients
noutparams    = Y.shape[1] # Number of parameters to predict
nsamples = X.shape[0]

if (perccoef < 0.99999):
    if nfrequencies%2==0:
        aks=(nfrequencies+2)/2
        bks=aks-2
    else:
        aks=(nfrequencies+1)/2
        bks=aks-1
    aks = int(aks)
    selectedaks = int(np.floor(aks * perccoef))
    indicestoselect = list(range(0,selectedaks)) + list(range(aks,aks+selectedaks-1))
    X = np.take(X,indicestoselect,axis=2)
    nfrequencies = X.shape[2]


# def normalizeinputmat(mat,newmins,newmaxs):

#     shapes = mat.shape
    
#     coefmins = np.full(newmins.shape,0.0)
#     coefmaxs = np.full(newmaxs.shape,0.0)

#     for coefind in range(0,shapes[1]):
#         for freq in range(0,shapes[2]):
#             coef = mat[:,coefind,freq]
#             coefmin = coef.min()
#             coefmax = coef.max()
#             coefmins[coefind,freq] = coefmin
#             coefmaxs[coefind,freq] = coefmax
#             if (abs(coefmax - coefmin) > 1e-15):#In case of normalization by column (bk removed in MATLAB)
#                 mat[:,coefind,freq] = newmins[coefind,freq] + (newmaxs[coefind,freq] - newmins[coefind,freq]) * (coef - coefmin) / (coefmax - coefmin)

#     return mat, coefmins, coefmaxs

# # different normalization: we normalize wrt to the min and max values of each
# # parameter at different samples
# def normalizeoutputmat(mat,newmins,newmaxs):
#     shapes = mat.shape
#     parammins = np.full(newmins.shape,0.0)
#     parammaxs = np.full(newmins.shape,0.0)
    
#     for i in range(0,shapes[1]):
#         param = mat[:,i]
#         parammin = param.min()
#         parammax = param.max()
#         parammins[i] = parammin
#         parammaxs[i] = parammax
#         mat[:,i] = newmins[i] + (newmaxs[i] - newmins[i]) * (param - parammin) / (parammax - parammin)

#     return mat, parammins, parammaxs



newmins = np.full([ncoefficients,nfrequencies],0.0)
newmaxs = np.full([ncoefficients,nfrequencies],1.0)

X,coefmins,coefmaxs = udl.normalizeinputmatDL(X,newmins,newmaxs,coefmins,coefmaxs)

newmins = np.full([noutparams],0.0)
newmaxs = np.full([noutparams],1.0)
Y,parammins,parammaxs = udl.normalizeoutputmatDL(Y,newmins,newmaxs,parammins,parammaxs)

Xtest = X
Ytest = Y

Ypred = dnn_model.predict(Xtest)

Xtest, _, _ = udl.normalizeinputmatDL(Xtest, coefmins, coefmaxs)
Ypred, _, _ = udl.normalizeoutputmatDL(Ypred, parammins, parammaxs) #Obtain real values, non normalized
Ytest, _, _ = udl.normalizeoutputmatDL(Ytest, parammins, parammaxs)

fig, axs = plt.subplots(2, 2)
axs[0,0].scatter(Ytest[:,0],Ypred[:,0])
axs[0,0].set_title('Left Ventricle Emax0')
axs[0,0].set_xlabel('real parameter')
axs[0,0].set_ylabel('predicted parameter')
axs[0,1].scatter(Ytest[:,1],Ypred[:,1])
axs[0,1].set_title('Left Ventricle EmaxRef0')
axs[0,1].set_xlabel('real parameter')
axs[0,1].set_ylabel('predicted parameter')
axs[1,0].scatter(Ytest[:,2],Ypred[:,2])
axs[1,0].set_title('Left Ventricle AGain_Emax')
axs[1,0].set_xlabel('real parameter')
axs[1,0].set_ylabel('predicted parameter')
axs[1,1].scatter(Ytest[:,3],Ypred[:,3])
axs[1,1].set_title('Left Ventricle kE')
axs[1,1].set_xlabel('real parameter')
axs[1,1].set_ylabel('predicted parameter')
axs[1,1].set_xlim(parammins[3],parammaxs[3])
axs[1,1].set_ylim(parammins[3],parammaxs[3])
fig.set_figheight(10)
fig.set_figwidth(15)
plt.savefig('DNN_Performance.eps')

#%%
filepath="/Users/jean.bonnemain/Documents/Code/0d_model/Modelica_Code/0D_Original/"
today = datetime.datetime.now()
outputfolder=os.getcwd() + '/' + today.strftime("%Y") + '_' + today.strftime("%m") + '_' + today.strftime("%d")

uo.prepareOutputFolder(outputfolder)

testsamples = Ypred.shape[0]
scipy.io.savemat('Ypred.mat',mdict={'Ypred': Ypred})
scipy.io.savemat('Ytest.mat',mdict={'Ytest': Ypred})

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