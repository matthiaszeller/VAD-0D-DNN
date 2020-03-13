from OMPython import OMCSessionZMQ
from random import seed
from random import random
import os
omc = OMCSessionZMQ()

from timeit import default_timer as timer
import utils_openmodelica as uo

import datetime
import numpy as np

filepath="/Users/jean.bonnemain/Documents/Code/0d_model/Modelica_Code/0D_Original/"
today = datetime.datetime.now()

outputfolder=os.getcwd() + '/' + today.strftime("%Y") + '_' + today.strftime("%m") + '_' + today.strftime("%d")

class ParameterPhysiological: 
    def __init__(self, name, minparam = None, maxparam = None, coef = None):
        self.minparam = minparam
        self.maxparam = maxparam
        self.name = name
        self.value = 0.
        self.randomsampling = True
        self.coef = coef
        
    def setValue(self, value):
        self.value = value
        self.randomsampling = False

    def sample(self):
        if self.randomsampling:
            self.value = self.minparam + (self.maxparam - self.minparam) * self.coef

for coef in np.arange(0,101,20):
    print(coef)
    param1 = ParameterPhysiological("Param_LeftVentricle_Emax0", 0.2, 2.95, coef/100)
    param2 = ParameterPhysiological("Param_LeftVentricle_EmaxRef0", 0.2, 2.392, coef/100)
    param3 = ParameterPhysiological("Param_LeftVentricle_AGain_Emax", 0.2, 0.475, coef/100)
    param4 = ParameterPhysiological("Param_LeftVentricle_kE", 0.011, 0.014, coef/100)
    listParameters = [param1, param2, param3, param4]
    numberofsamples = 1
    uo.prepareOutputFolder(outputfolder)
    suffix = str(coef)
    uo.launchSimulation(filepath, listParameters, suffix, outputfolder)

