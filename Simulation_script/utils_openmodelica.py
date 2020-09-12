from OMPython import OMCSessionZMQ
from random import seed
from random import random
import os
from setup import *
import numpy.random

from timeit import default_timer as timer


class Parameter: 
    def __init__(self, name, minparam = None, maxparam = None, epsilon=None):
        self.minparam = minparam
        self.maxparam = maxparam
        self.name = name
        self.value = 0
        self.randomsampling = True
        self.epsilon = epsilon
        
    def setValue(self, value):
        self.value = value
        self.randomsampling = False

    def sample(self):
        if self.randomsampling:
            self.value = self.minparam + (self.maxparam - self.minparam) * random()

    def sample_epsilon(self):
        if self.randomsampling:
            self.value = ((self.maxparam + self.minparam)/2) * (1+self.epsilon * numpy.random.normal(0,1))


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


def sampleParamsEpsilon(param_lst):
    for p in param_lst:
        p.sample_epsilon()


def runSimulation(N, param_lst, output_folder, file, LVAD, samplingfun=sampleParams,
                  om_sim_settings=None, override_params=None):
    """Run the whole simulation. The process involves two steps:
    1- Build the model. It creates an executable in the build folder.
    2- Run N simulations (does not require compilation).
    For each simulation, generate new parameter values and override them in the original model,
    according to the parameters specified in param_lst.
    Note: we can run simulation without build by using ./exec_name -override=paramName=ParamValue

    :param dict[str, any] override_params: dictionnary of 0D model parameters to override
    :param dict[str, any] om_sim_settings: dictionnary of OpenModelica simulation settings
    """
    # 1. Create folders
    try:
        out = prepareOutputFolder(output_folder)
    except Exception:
        raise Exception("Could not create the folder '{}' for unknown reason".format(output_folder))
    # Do not overwrite in existing folder ! Abort the program
    if out == False:
        raise ValueError("The folder '{}' already exists".format(output_folder))
    # 2. Run OMC session
    omc = OMCSessionZMQ()
    # 3. Determine model name
    model_name = "Mathcard.Applications.Ursino1998.Ursino1998Model"
    if LVAD: model_name = "Mathcard.Applications.Ursino1998.HMIII.Ursino1998Model_VAD2"

    # 4. Build the model (done only once)
    # An executable will be created,
    env = {}
    runcmd(omc, "loadModel(Modelica)", env)
    print("loadFile(\"{}\")".format(file))
    runcmd(omc, "loadFile(\"{}\")".format(file), env)
    runcmd(omc, "instantiateModel({})".format(model_name), env)

    cmd_compile = "simulate({}, stopTime=30.0, numberOfIntervals=2000, " \
                  "simflags=\"-emit_protected\", outputFormat=\"mat\")".format(model_name)
    if om_sim_settings is not None:
        cmd_compile = f'simulate({model_name}, ' + ', '.join([
            k + '=' + str(v)
            for k,v in om_sim_settings.items()
        ]) + ')'

    runcmd(omc, cmd_compile, env)

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

    override_params = ','.join(k + '=' + str(v) for k, v in override_params.items()) \
                      if override_params is not None else ''

    # Loop
    for n in range(N):
        # Sample the parameters
        samplingfun(param_lst)
        # Save parameter values
        dic = {p.name : p.value for p in param_lst}
        param_data[n] = dic
        # Create the command to override parameters
        override_cmd = override_cmd_template.format(*[p.value for p in param_lst])
        # Specify output file
        output_file = output_file_template.format(n)
        # Simulate without build
        cmd = "./" + model_name + " -override="+override_cmd+override_params \
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


def runTestSimulation(Ytest, Ytest_pred, param_lst, output_dnn_test, modelica_file_path):
    """Given the predicted responses of the neural network and the exact ones (those used to generate input data for DNN)
    run simulations to compare hemodynamics between the two."""
    # ======== PREPARE
    try:
        out = prepareOutputFolder(output_dnn_test)
    except Exception:
        print("\n===================== FATAL ERROR =====================")
        print("Could not create the folder '{}' for unknown reason".format(output_dnn_test))
    # Do not overwrite in existing folder ! Abort the program
    if out == False:
        print("\n===================== FATAL ERROR =====================")
        print("The folder '{}' already exists".format(output_dnn_test))
        exit()

    # ======== PREPARE SESSION
    omc = OMCSessionZMQ()
    model_name = "Mathcard.Applications.Ursino1998.Ursino1998Model"
    if SIMULATION_LVAD: model_name = "Mathcard.Applications.Ursino1998.HMIII.Ursino1998Model_VAD2"

    # ======== LOAD FILES
    env = {}
    runcmd(omc, "loadModel(Modelica)", env)
    runcmd(omc, "loadFile(\"{}\")".format(modelica_file_path), env)
    runcmd(omc, "instantiateModel({})".format(model_name), env)
    runcmd(omc, "simulate({}, stopTime=30.0, numberOfIntervals=2000, \
    simflags=\"-emit_protected\", outputFormat=\"mat\")".format(model_name), env)

    # ======== PREPARE SIMULATION
    override_cmd_template = [p.name+"={}" for p in param_lst]
    override_cmd_template = ','.join(override_cmd_template)
    # Prepare output file format
    # Suffix (after simulation number) will be 'predicted' or 'exact'
    output_file_template = out + model_name.split('.')[-1] + "_output_{}_{}.mat"

    # ========= SIMULATION LOOP
    for i in range(Ytest.shape[0]):
        exact_params = Ytest[i, :]
        pred_params = Ytest_pred[i, :]

        # == Launch simulation with exact parameters
        override_cmd = override_cmd_template.format(*exact_params)
        output_file = output_file_template.format(i, 'exact')
        # Simulate without build
        cmd = "./" + model_name + " -override="+override_cmd \
                   + " -r="+output_file + " -emit_protected"
        q("<RUNNING> " + cmd)
        start = timer()
        res = os.system(cmd)
        end = timer()
        q("Simulation time (n={}, exact): {:.5}".format(i, end-start))

        # == Launch simulation with predicted parameters
        override_cmd = override_cmd_template.format(*pred_params)
        output_file = output_file_template.format(i, 'predicted')
        # Simulate without build
        cmd = "./" + model_name + " -override="+override_cmd \
                   + " -r="+output_file + " -emit_protected"
        q("<RUNNING> " + cmd)
        start = timer()
        res = os.system(cmd)
        end = timer()
        q("Simulation time (n={}, precited): {:.5}".format(i, end-start))