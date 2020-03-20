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
    shutil.move(filepath + outputfile, outputfolder + "/models/" + outputfile)

    
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


def launchSimulation(filepath, listParameters, suffix, outputfolder, LVAD=False):
    # Store the output of the commands
    env = {}

    def run(cmd):
        q("<RUNNING> " + cmd)
        env[cmd] = omc.sendExpression(cmd)
        # Print the output, except for instanciateModel (too verbose)
        if not "instantiateModel" in cmd:
            q(str(env[cmd]))

    # 1. DETERMINE OUTPUT FILE
    outputfile = "Mathcard" + "_" + suffix + ".mo"

    # 2. CREATE NEW MODELICA FILE AND CHANGE PARAMS
    changeValueInFile(filepath, listParameters, outputfile, outputfolder)

    # 3. LOAD MODELS
    model_name = "Mathcard.Applications.Ursino1998.Ursino1998Model"
    if LVAD: model_name = "Mathcard.Applications.Ursino1998.HMIII.Ursino1998Model_VAD2"

    run("loadModel(Modelica)")
    run("loadFile(\"" + outputfolder + "/models/" + outputfile + "\")")
    run("instantiateModel({})".format(model_name))

    # 4. RUN THE SIMULATION
    start = timer()
    run("simulate({}, \
    stopTime=20.0, numberOfIntervals=500, \
    simflags=\"-emit_protected\", outputFormat=\"mat\")".format(model_name))
    end = timer()

    # 5. REPORT SIMULATION TIME
    print("Simulation time [{}]\t{:.4}".format(suffix, end - start))

    # 6. MOVE RESULT FILE
    output_name = "Ursino1998Model"
    if LVAD: output_name += "_VAD2"
    matrixoutput = output_name + "_output_" + suffix + ".mat"
    q("shutil.move: {}_res.mat -> ".format(model_name) + outputfolder + "/outputs/" + matrixoutput)
    shutil.move(model_name + "_res.mat", outputfolder + "/outputs/" + matrixoutput)

