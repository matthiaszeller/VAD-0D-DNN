from OMPython import OMCSessionZMQ
from random import seed
from random import random
import os
omc = OMCSessionZMQ()

from timeit import default_timer as timer
import utils_openmodelica as uo

import datetime

# File paths
from setup import *

today = datetime.datetime.now()

#outputfolder=os.getcwd() + '/' + today.strftime("%Y") + '_' + today.strftime("%m") + '_' + today.strftime("%d")
outputfolder += "_" + today.strftime("%Y") + '_' + today.strftime("%m") + '_' + today.strftime("%d")

q("<INFO> The output folder is "+outputfolder)

# list of parameters which we vary
param1 = uo.Parameter("Param_LeftVentricle_Emax0", 0.2, 2.95)
param2 = uo.Parameter("Param_LeftVentricle_EmaxRef0", 0.2, 2.392)
param3 = uo.Parameter("Param_LeftVentricle_AGain_Emax", 0.2, 0.475)
param4 = uo.Parameter("Param_LeftVentricle_kE", 0.011, 0.014)

listParameters = [param1, param2, param3, param4]

numberofsamples = 1

outputCreated = uo.prepareOutputFolder(outputfolder)
#if not outputCreated:
#    print("<INFO> The program is aborted.")
#    exit(0)

for indexsample in range(0,numberofsamples):
    suffix = str(indexsample)
    uo.launchSimulation(filepath, listParameters, suffix, outputfolder)
