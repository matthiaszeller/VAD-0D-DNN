from OMPython import OMCSessionZMQ
from random import seed
from random import random
import os
from setup import *
import shutil

from timeit import default_timer as timer

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
    if os.path.exists(outputfolder):
        return False

    # Create main folder
    os.mkdir(outputfolder)
    # Create folder for output data
    data_output = outputfolder + "/outputs/"
    os.mkdir(data_output)
    q("<INIT> Output folder: " + outputfolder)
    return data_output

def runcmd(omc, cmd, env):
    """Run openmodelica command with the OMCSessionZMQ instance passed as an argument"""
    # Debug: display command to be run
    q("<RUNNING> " + cmd)
    # Run the actual command and save result
    env[cmd] = omc.sendExpression(cmd)
    # Print (debug mode) the output, except for instanciateModel (too verbose)
    if not "instantiateModel" in cmd:
        q("<CMD OUTPUT> " + str(env[cmd]))


def sampleParams(param_lst):
    for p in param_lst:
        p.sample()


def runSimulation(N, param_lst, output_folder, file, LVAD):
    """Run the whole simulation. The process involves two steps:
    1- Build the model. It creates an executable in the build folder.
    2- Run N simulations (does not require compilation).
    For each simulation, generate new parameter values and override them in the original model,
    according to the parameters specified in param_lst.
    Note: we can run simulation without build by using ./exec_name -override=paramName=ParamValue"""
    # 1. Create folders
    try:
        out = prepareOutputFolder(output_folder)
    except Exception:
        print("\n===================== FATAL ERROR =====================")
        print("Could not create the folder '{}' for unknown reason".format(output_folder))
    # Do not overwrite in existing folder ! Abort the program
    if out == False:
        print("\n===================== FATAL ERROR =====================")
        print("The folder '{}' already exists".format(output_folder))
        exit()
    # 2. Run OMC session
    omc = OMCSessionZMQ()
    # 3. Determine model name
    model_name = "Mathcard.Applications.Ursino1998.Ursino1998Model"
    if LVAD: model_name = "Mathcard.Applications.Ursino1998.HMIII.Ursino1998Model_VAD2"
    ###############################
    # 4. Build the model (done only once)
    # An executable will be created,
    env = {}
    runcmd(omc, "loadModel(Modelica)", env)
    print("loadFile(\"{}\")".format(file))
    runcmd(omc, "loadFile(\"{}\")".format(file), env)
    runcmd(omc, "instantiateModel({})".format(model_name), env)
    runcmd(omc, "simulate({}, stopTime=20.0, numberOfIntervals=500, \
    simflags=\"-emit_protected\", outputFormat=\"mat\")".format(model_name), env)
    ###############################
    # 5. Prepare the simulation
    # Prepare parameter recording
    param_data = {}
    # Prepare parameter overriding
    # Idea: use override_cmd_template with format to easily replace the param values
    override_cmd_template = [p.name+"={}" for p in param_lst]
    override_cmd_template = ','.join(override_cmd_template)
    # Prepare output file format
    # Note: we don't take the whole model name, only the end of it (-> use split)
    output_file_template = out + model_name.split('.')[-1] + "_output_{}.mat"
    ###############################
    # Loop
    for n in range(N):
        # Sample the parameters
        sampleParams(param_lst)
        # Save parameter values
        dic = {p.name : p.value for p in param_lst}
        param_data[n] = dic
        # Create the command to override parameters
        override_cmd = override_cmd_template.format(*[p.value for p in param_lst])
        # Specify output file
        output_file = output_file_template.format(n)
        # Simulate without build
        cmd = "./" + model_name + " -override="+override_cmd \
                   + " -r="+output_file + " -emit_protected"
        q("<RUNNING> " + cmd)
        start = timer()
        res = os.system(cmd)
        end = timer()
        q("Simulation time (n={}): {:.5}".format(n, end-start))
    writeParamData(param_data, output_folder)

def writeParamData(data, output_folder):
    file_path = output_folder + "/" + "parameters.txt"
    # Format the file content in a string
    cols = [name for name in data[0]]
    header = ','.join(['n'] + cols)
    txt = header
    for n, values in data.items():
        line = str(n)
        for col in cols: line += ","+str(values[col])
        txt += "\n" + line
	# Open the file and write
    with open(file_path, 'w') as f:
        f.write(txt)
