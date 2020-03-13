from OMPython import OMCSessionZMQ
from random import seed
from random import random
import os
from setup import *
import shutil

from timeit import default_timer as timer

omc = OMCSessionZMQ()


class Parameter: 
    def __init__(self, name, minparam = None, maxparam = None):
        self.minparam = minparam
        self.maxparam = maxparam
        self.name = name
        self.value = 0
        self.randomsampling = True
        
    def setValue(self, value):
        self.value = value
        self.randomsampling = False

    def sample(self):
        if self.randomsampling:
            self.value = self.minparam + (self.maxparam - self.minparam) * random()

def changeValueInFile(filepath, listParameter, outputfile, outputfolder):
    filename = filepath + "Mathcard"
    numparams = len(listParameter)
    numParameters = -1
    count = 1
    filehandler = open(filename + ".mo", "r");
    filehandlerout = open(filepath + "/" + outputfile, "w")
    for line in filehandler:
        copyLine = 1
        # Principle: find the following model in the file,
        # once it is found, we modify as many as `numparams` parameters,
        # thus it's ensured that only the params of that model are modified
        if "model ModelParametersNH" in line:
            numParameters = 0
        for parameter in listParameter:
            if parameter.name in line and numParameters < numparams and numParameters >= 0:
                indexequal = line.rfind("=")
                line = line[:indexequal]
                parameter.sample()
                line += "= " + str(parameter.value) + ";\n"
                filehandlerout.write(line)
                numParameters = numParameters + 1
                copyLine = 0
        if copyLine:
            filehandlerout.write(line)
        count = count + 1
    filehandler.close()
    filehandlerout.close()

    q("shutil.move :" + filepath + outputfile + "->" + outputfolder + "/models/" + outputfile)
    #os.rename(filepath + "/" + outputfile, outputfolder + "/models/" + outputfile)
    shutil.move(filepath + "/" + outputfile, outputfolder + "/models/" + outputfile)
    #shutil.copy(filepath + "/" + outputfile, outputfolder + "/models/" + outputfile)
    #shutil.rmtree(filepath + "/" + outputfile)

    
def prepareOutputFolder(outputfolder):
    try:
        q("1")
        os.mkdir(outputfolder)
        q("2")
        os.mkdir(outputfolder + '/models')
        q("3")
        os.mkdir(outputfolder + '/outputs')
        q("prepareOutputFolder: models and outputs folders were created in"+ outputFolder)
        return True
    except Exception:
        print("<ERROR> Directory " + outputfolder + " exists!")
    return False
    
def launchSimulation(filepath, listParameters, suffix, outputfolder):
    outputfile = "Mathcard" + "_" + suffix + ".mo"
    changeValueInFile(filepath, listParameters, outputfile, outputfolder)
    omc.sendExpression("loadModel(Modelica)")
    omc.sendExpression("loadFile(\"" + outputfolder + "/models/" + outputfile + "\")")
    omc.sendExpression("instantiateModel(Mathcard.Applications.Ursino1998.Ursino1998Model)")
    matrixoutput = "Ursino1998Model_output" "_" + suffix + ".mat"
    start = timer()
    omc.sendExpression("simulate(Mathcard.Applications.Ursino1998.Ursino1998Model, stopTime=20.0, numberOfIntervals=500, simflags=\"-emit_protected\")")
    end = timer()
    print("Simulation time [{}]\t{:.4}".format(suffix, end - start))
    shutil.move("Mathcard.Applications.Ursino1998.Ursino1998Model_res.mat", outputfolder + "/outputs/" + matrixoutput)

