
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

# ============ IMPLEMENTATION =============

# One can run the script with arguments to run specific parts
run_train, run_test, path = udl.manage_args(sys.argv)

print('Tensorflow version:', tf.__version__)

if run_train:
    # PERCENTAGE OF COEFFICIENTS TO KEEP
    perccoef = 0.05
    # Launch training
    udl.train_dnn(perccoef, files_path=path)

print("\nExiting to keep control...")
exit()


# cut the input matrix (e.g. to keep only a subset of the frequencies)
# freqmax = 100
# X = X[:,:,0:freqmax]

# for simplicity this is specific to a three dimensional tensor in which the
# the second component differentiates between aks and bks
def normalizeinputmat(mat,newmins,newmaxs):
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

    return mat, coefmins, coefmaxs

# different normalization: we normalize wrt to the min and max values of each
# parameter at different samples
def normalizeoutputmat(mat,newmins,newmaxs):
    shapes = mat.shape
    parammins = np.full(newmins.shape,0.0)
    parammaxs = np.full(newmins.shape,0.0)
    
    for i in range(0,shapes[1]):
        param = mat[:,i]
        parammin = param.min()
        parammax = param.max()
        parammins[i] = parammin
        parammaxs[i] = parammax
        mat[:,i] = newmins[i] + (newmaxs[i] - newmins[i]) * (param - parammin) / (parammax - parammin)

    return mat, parammins, parammaxs


newmins = np.full([ncoefficients,nfrequencies],0.0)
newmaxs = np.full([ncoefficients,nfrequencies],1.0)

X,coefmins,coefmaxs = udl.normalizeinputmatDL(X,newmins,newmaxs)
np.save("coefmins",coefmins)
np.save("coefmaxs",coefmaxs)

newmins = np.full([noutparams],0.0)
newmaxs = np.full([noutparams],1.0)
Y,parammins,parammaxs = udl.normalizeoutputmatDL(Y,newmins,newmaxs)
np.save("parammins",parammins)
np.save("parammaxs",parammaxs)


# test and validation
perctest = 0.05
perctrain = 1 - perctest
percvalidation = 0.2
samplestrain = int(np.floor(perctrain * nsamples))
samplestest = nsamples - samplestrain

Xtrain = X[0:samplestrain,:,:]
Ytrain = Y[0:samplestrain,:]
Xtest = X[samplestrain+1:,:,:]
Ytest = Y[samplestrain+1:,:]

model = keras.Sequential([
    keras.layers.Flatten(input_shape=(ncoefficients, nfrequencies)),
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

history = model.fit(Xtrain, Ytrain, epochs=1000, validation_split=percvalidation)

plt.semilogy(history.history['loss'], label = 'loss')
plt.semilogy(history.history['val_loss'], label = 'val_loss')
plt.legend()
plt.savefig('Losses.eps')

#Save the model in the .h5 format
model.save('DNN_0D_Model.h5')

Ypred = model.predict(Xtest)

Xtest, _, _ = normalizeinputmat(Xtest, coefmins, coefmaxs)
Ypred, _, _ = normalizeoutputmat(Ypred, parammins, parammaxs) #Obtain real values, non normalized
Ytest, _, _ = normalizeoutputmat(Ytest, parammins, parammaxs)
print(Ypred[0,:])
print(Xtest[0,0,:]) #return frequency for systemic arteries
#scipy.io.savemat('Ypred.mat',mdict={'Ypred': Ypred})
#scipy.io.savemat('Ypred.mat',mdict={'Xtest': Xtest})


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