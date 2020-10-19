"""Use this to run simulations with exact vs predicted parameters,
once a DNN architecture has been found to accurately predict heart failure parameters.

You should have launched `script_pipeline.py` first.

Run this script from the folder `dnns`, whose contents should look like:
    .
    ├── dnn_10_layers_128_neurons
    ├── dnn_10_layers_16_neurons
    ├── [...]
    ├── dnn_8_layers_32_neurons
    ├── dnn_8_layers_64_neurons
    └── normdata.mat

"""

from script_pipeline import om_build_settings, om_runtime_args, \
    project_path, modelica_file_path, gen_trash_folder, list_parameters
from OMPython import OMCSessionZMQ
import sys
import os
from os import path
from time import time
import argparse
from scipy.io import loadmat
import numpy as np
import shutil

sys.path.insert(1, path.join(project_path, 'Simulation_script/'))
import utils_openmodelica as uo


SIMULATION_LVAD = True


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('dnn_folder', help='folder containing trained DNN data', type=str)
    args = parser.parse_args()

    if args.dnn_folder[-1] == '/':
        args.dnn_folder = args.dnn_folder[:-1]

    print('Loading normalized data to extract rpms:')
    norm = loadmat('normdata.mat')
    rpms = norm['Xtest'][:, -1]
    # Unnormalize
    rpms = rpms * (norm['Xmaxs'][:, -1] - norm['Xmins'][:, -1]) + norm['Xmins'][:, -1]
    del norm
    print(f'Rpms extracted (shape {rpms.shape}):\n{rpms}')

    print('Loading Ytest.txt and Ytestpred.txt from', args.dnn_folder)
    Ytest = np.genfromtxt(os.path.join(args.dnn_folder, 'Ytest.txt'))
    Ytesthat = np.genfromtxt(os.path.join(args.dnn_folder, 'Ytestpred.txt'))

    output_folder = f'../{args.dnn_folder}_test'

    print('Generating trash folder...')
    trash = gen_trash_folder('.')
    os.chdir(trash)

    print(f'Running simulations, PID = {os.getpid()}...')
    try:
        run_test_simulations(Ytest, Ytesthat, rpms, list_parameters, output_folder,
                             modelica_file_path, om_build_settings, om_runtime_args)
    except Exception as e:
        print(e)
    finally:
        print('Deleting trash folder...')
        os.chdir('..')
        shutil.rmtree(trash)


def run_test_simulations(Ytest, Ytest_pred, rpms, param_lst, output_dnn_test, modelica_file_path,
                         om_build_settings, om_runtime_args):
    """Given the predicted responses of the neural network and the exact ones (those used to generate input data for DNN)
    run simulations to compare hemodynamics between the two."""
    # ======== PREPARE
    try:
        out = uo.prepareOutputFolder(output_dnn_test)
    except Exception:
        raise Exception("Could not create the folder '{}' for unknown reason".format(output_dnn_test))
    # Do not overwrite in existing folder ! Abort the program
    if out == False:
        raise ValueError("The folder '{}' already exists".format(output_dnn_test))

    # ======== PREPARE SESSION
    omc = OMCSessionZMQ()
    model_name = "Mathcard.Applications.Ursino1998.Ursino1998Model"
    if SIMULATION_LVAD: model_name = "Mathcard.Applications.Ursino1998.HMIII.Ursino1998Model_VAD2"

    # ======== LOAD FILES
    env = {}
    log = print
    uo.runcmd(omc, "loadModel(Modelica)", env, log)
    uo.runcmd(omc, "loadFile(\"{}\")".format(modelica_file_path), env, log)
    uo.runcmd(omc, "instantiateModel({})".format(model_name), env, log)

    cmd_compile = "simulate({}, stopTime=30.0, numberOfIntervals=2000, " \
                  "simflags=\"-emit_protected\", outputFormat=\"mat\")".format(model_name)
    if om_build_settings is not None:
        cmd_compile = f'simulate({model_name}, ' + ', '.join([
            k + '=' + str(v)
            for k,v in om_build_settings.items()
        ]) + ')'
    uo.runcmd(omc, cmd_compile, env, log)

    # ======== PREPARE SIMULATION
    override_cmd_template = [p.name + "={}" for p in param_lst] + ['Param_LVAD_RPM={}']
    override_cmd_template = ','.join(override_cmd_template)
    # Prepare output file format
    # Suffix (after simulation number) will be 'predicted' or 'exact'
    output_file_template = out + model_name.split('.')[-1] + "_output_{}_{}.mat"

    # ========= SIMULATION LOOP
    for i in range(Ytest.shape[0]):
        exact_params = Ytest[i, :]
        pred_params = Ytest_pred[i, :]
        rpm = rpms[i]

        # == Launch simulation with exact parameters
        override_cmd = override_cmd_template.format(*exact_params, rpm)
        output_file = output_file_template.format(i, 'exact')
        # Simulate without build
        cmd = "./" + model_name + " -override=" + override_cmd \
              + " -r=" + output_file + " -emit_protected " + om_runtime_args
        log("<RUNNING> " + cmd)
        start = time()
        res = os.system(cmd)
        end = time()
        log("Simulation time (n={}, exact): {:.5}".format(i, end - start))

        # == Launch simulation with predicted parameters
        override_cmd = override_cmd_template.format(*pred_params, rpm)
        output_file = output_file_template.format(i, 'predicted')
        # Simulate without build
        cmd = "./" + model_name + " -override=" + override_cmd \
              + " -r=" + output_file + " -emit_protected " + om_runtime_args
        log("<RUNNING> " + cmd)
        start = time()
        res = os.system(cmd)
        end = time()
        log("Simulation time (n={}, precited): {:.5}".format(i, end - start))


if __name__ == '__main__':
    main()
