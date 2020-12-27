"""This script wraps all the code related to:
    - simulation setup
    - simulation data generation
    - data pre-processing
    - DNN training

When running the script, you are first prompted to continue.
Advice:
    - to do some tests, set `N_samples` to a very low value (but this might trigger error for DNN training)
    - read the configuration displayed when you run the script
    - use `script_stats_dnn.py` to find the best architectures before running this script
    - you can use `script_check_sim_data.py` to check consistency of simulation data

General script notions:
    - A dataset is made of `N_samples` simulation files
    - A dataset is characterized by its global simulation settings: pump speed (RPM),
      whether the LVAD is used, and whether artificial pulse (AP) is activated
    - The word "pipeline" used in this script refers to generating data for a single dataset
    - This script runs `N_processes` processes in parallel and assigns 1 pipeline to each process
    - You can run the script several times with the same `root_output_folder`,
      but don't change `N_samples` in that case
    - If you run this script in a folder with existing data, it tries to resume the pipeline
      and outputs an error if the pipeline cannot be resumed
    - Do NOT rename/move folders nor any file created by the script

Example of usage:
    1. set `configurations` with only 1 DNN architecture (you can train other DNNs later),
       and at most `N_processes` configurations
    2. change other script variables according to what you need (in `SETUP` section)
    3. run this script
    4. repeat from 1. and set the other configurations you want (don't change N_samples and provide same number of
       coefficients)
    5. once all simulation data is computed, you can now train additional DNN configurations,
       set the `configurations` variable accordingly
    5. re-run the script - it trains the additional DNNs (0D data is not computed again)

Warning:
    * You must edit `setup_preprocessing.py` according to settings you chose
      in this script (time discretization according to values in `om_build_settings`)

Once you identified a suitable DNN architecture, you can run script_test_dnn.py.
"""

# TODO: a huge bottleneck seems to be pre-processing when using an HDD:
#       very (!) inefficient if all processes read the N_samples at same time


# ======================================================= #
# ----------------------- IMPORTS ----------------------- #
# ======================================================= #

from os import path
import os
import shutil
import multiprocessing
import random
import json
import subprocess
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import argparse
import pickle
import time
from datetime import datetime
import scipy.io as sio

# Reduce Tensorflow verbosity
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '2'

# --- Import project modules
current_script_path = path.dirname(path.abspath(__file__))
project_path = path.dirname(current_script_path)

sys.path.insert(1, path.join(project_path, 'Simulation_script/'))
import utils_openmodelica as uo

sys.path.insert(1, path.join(project_path, 'pipelining/'))
from om_server_handler import OMServerHanlder
from Logger import Logger

sys.path.insert(1, path.join(project_path, 'Deep_learning/'))
import utils_deeplearning as ud

# ========================================================= #
# ------------------------- SETUP ------------------------- #
# ========================================================= #

# -------------------- GENERAL SETUP

# --- Root folder containing all data
root_output_folder = '/media/maousi/Raw/lvad'
root_trash_folder = path.join(root_output_folder, 'trash')

# --- Multiprocessing
N_processes = multiprocessing.cpu_count()

# --- Pump configurations (pump speed, LVAD, artificial Pulse, DNN architectures)
# DNN architectures: list of tuples (hidden_layers, neurons, aks)
configurations = [
    (rpm, True, True)
    for rpm in range(4600, 6100, 100)
]

# --- DNN architectures (hidden layers, neurons)
architectures = [
    (layers, neurons)
    for layers in [8, 10, 12]
    for neurons in [16, 32, 64, 128]
] + [
    (layers, neurons)
    for neurons in [32, 64, 128, 256]
    for layers in [3, 4, 5, 6]
]

# -------------------- SIMULATION SETUP

# The number of simulations in a single dataset
N_samples = 10000

# OpenModelica build settings: simulation time, number of intervals, etc...
om_build_settings = {
    'stopTime': 30.0,
    'numberOfIntervals': 2000,
    'outputFormat': '"mat"',
    'simflags': '"-emit_protected"',
}

# 0D Model parameters overriding
global_0D_params_override = {

}

# Arguments with which to run executables
om_runtime_args = '-lv=-LOG_SUCCESS' # reduce verbosity

# The path to `Mathcard.mo`
modelica_file_path = path.join(project_path, 'modelica/original/Mathcard.mo')

# -------------------- 0D-model parameters that vary during simulation
param1 = uo.Parameter("Param_LeftVentricle_Emax0", 0.2, 2.95)
param2 = uo.Parameter("Param_LeftVentricle_EmaxRef0", 0.2, 2.392)
param3 = uo.Parameter("Param_LeftVentricle_AGain_Emax", 0.2, 0.475)
param4 = uo.Parameter("Param_LeftVentricle_kE", 0.011, 0.014)
list_parameters = [param1, param2, param3, param4]


# -------------------- SETUP CHECKUP
if path.isfile(root_output_folder):
    raise ValueError('Invalid root_output_folder')
elif not path.exists(root_output_folder):
    os.mkdir(root_output_folder)

# ---------- PREPROCESSING
preprocessing_script_path = path.join(project_path, 'Pre-processing', 'script_createdataset.py')


# -------------------- MINIMAL REQUIRED DNN FOLDER CONTENTS
dnn_folder_contents = [ # TODO: revise minimal required folder contents
    'history.bin', 'DNN_0D_Model.h5', 'Xtest_norm.npy',
    'Ytest_norm.npy', 'parammins.npy', 'parammaxs.npy'
]

# ========================================================= #
# ---------------------- SCRIPT BODY ---------------------- #
# ========================================================= #


# --------------------------------------------------------- #
#                     DATASET SETTINGS                      #
# --------------------------------------------------------- #


def format_folder_name(configuration):
    rpm, lvad, artpulse = configuration
    lvad = '_LVAD' if lvad else ''
    artpulse = '_AP' if artpulse else ''
    return f'{rpm}{lvad}{artpulse}'


def parse_folder_name(name):
    """Returns RPM of dataset"""
    return int(name.split('_')[0])


def build_0D_parameters_from_config(configuration):
    rpm, lvad, artpulse = configuration
    dic = {}
    # --- Pump speed
    dic['Param_LVAD_RPM'] = rpm
    # --- Artificial Pulse
    if not artpulse:
        dic['HMIII_Pulse_Amplitude'] = 0

    return dic


def format_dnn_folder_name(architecture):
    layers, neurons = architecture
    return f'dnn_{layers}_layers_{neurons}_neurons'


# ------------------------------------------------------- #
#                        SCRIPTING                        #
# ------------------------------------------------------- #


def prompt_initial_validate():
    def display(var, dic, indent=30):
        n = '\t\t'
        if isinstance(var, str):
            n = indent - len(var)
        print(f' - {var}:{n * " "}{dic[var]}')

    sep = '=' * 20
    print('\n' + sep)
    print('PIPELINING SCRIPT')
    print('\nPlease, carefully verify the settings below before proceeding:')
    prompt_vars = ['N_processes', 'N_samples', 'root_output_folder',
                   'om_build_settings', 'global_0D_params_override']

    for var in prompt_vars:
        display(var, globals())

    print(' - configurations:')
    print('\tRPM,   LVAD, AP')
    for c in configurations:
        print('\t' + str(c))

    i = input('\nDo you want to continue ? [y/N] ')
    if i != 'y':
        exit()


def main():
    global root_output_folder

    # ------------- SCRIPT ARGUMENTS
    parser = argparse.ArgumentParser()
    parser.add_argument('-l', '--logs',
                        help='file name in which to concatenate logs, write in file and exit',
                        type=str, default=None)
    parser.add_argument('-p', '--path',
                        help='override the root folder',
                        type=str, default=None)
    parser.add_argument('-s', '--status',
                        help='print status of each dataset before running the script',
                        action='store_true')
    parser.add_argument('--concat', action='store_true',
                        help='skip step 1 & 2 (concat datasets and train DNNs)')
    parser.add_argument('--train', action='store_true',
                        help='skip step 1, 2, 3 (train DNNs)')
    parser.add_argument('aks', type=int,
                        help='number of aks coeffs to select (control DNN input size)')
    args = parser.parse_args()

    if args.path is not None:
        root_output_folder = args.path
    if args.logs:
        print('Concatenating logs...')
        df = concat_logs(root_output_folder)
        save_concat_logs(args.logs, df)
        print('Exiting...')
        return

    # ------------- USER INTERACTION
    if args.status:
        # Print status of each dataset
        class Dummy:
            def log(self, msg):
                print(msg)

        print('Pipeline status:')
        for folder in os.listdir(root_output_folder):
            folder_path = os.path.join(root_output_folder, folder)
            if 'trash' in folder or not os.path.isdir(folder_path):
                continue

            status = check_pipeline_status(folder_path, Dummy())
            print('\t', folder, status)

    prompt_initial_validate()

    # ------------- PREPARATION
    mainlogger = Logger(path.join(root_output_folder, 'mainlog.log'),
                        print=True, buffer_size=0)
    mainlogger.log('Initialization...')

    if not args.concat and not args.train:
        if not path.exists(root_trash_folder):
            mainlogger.log('Creating trash folder...')
            os.mkdir(root_trash_folder)
        else:
            mainlogger.log('Trash folder detected')

        # Run a thread in background to kill OpenModelica sessions that never close
        mainlogger.log('Starting OMServerHanlder...')
        server_handler = OMServerHanlder(logfun=mainlogger.log)

        # ------------- RUN PIPELINES
        mainlogger.log('Launching processes from pool of size ' + str(N_processes))
        ps_pool = multiprocessing.Pool(processes=N_processes)
        ps_pool.map(pipeline, configurations)

        server_handler.kill_sessions()
        # Remove trash files
        if path.exists(root_trash_folder):
            mainlogger.log('Deleting trash data...')
            shutil.rmtree(root_trash_folder)

    # ------------- Concatenate datasets
    os.chdir(root_output_folder)
    if not args.train:
        root_files = os.listdir('.')
        # Keep folders
        root_files = filter(lambda fp: os.path.isdir(fp), root_files)
        # Don't keep trash
        root_files = filter(lambda fp: 'trash' not in fp, root_files)
        # Keep only complete datasets
        root_files = filter(lambda fp: check_pipeline_status(fp, mainlogger) == 2, root_files)

        root_files = list(root_files)
        datasets = {parse_folder_name(folder): folder for folder in root_files}

        oldstdout = sys.stdout
        sys.stdout = mainlogger
        mainlogger.set_stdout(oldstdout)
        X, Y = ud.generate_dataset(datasets, args.aks, perc_coef=None)
        # Save
        dataset_description = {
            'RPMs': list(datasets.keys()),
            'aks': args.aks
        }
        sio.savemat('dataset.mat', mdict={
            'description': dataset_description,
            'X': X,
            'Y': Y,
        })
        sys.stdout = oldstdout
        mainlogger.set_stdout(None)

    # ------------- TRAIN DNNS
    # in this case, the script basically starts here -> load data
    mainlogger.log('Loading existing data "dataset.mat", description:')
    data = sio.loadmat('dataset.mat')
    X, Y = data['X'], data['Y']
    mainlogger.log(str(data['description']))
    gen_training_data(X, Y, mainlogger)

    oldstdout = sys.stdout
    sys.stdout = mainlogger
    mainlogger.set_stdout(oldstdout)
    ps_pool = multiprocessing.Pool(processes=2)
    ps_pool.map(train_dnn, architectures)
    sys.stdout = oldstdout
    mainlogger.set_stdout(None)

    # Concatenate log files
    mainlogger.log('Concatenating log files...')
    df = concat_logs(root_output_folder)
    save_concat_logs(path.join(root_output_folder, 'logs.log'), df)

    mainlogger.flush()


# ------------------------------------------------------- #
#                      DATA PIPELINE                      #
# ------------------------------------------------------- #


def gen_trash_folder(root_trash_folder, configuration=None):
    """Create a folder with a unique name to perform temporary computations.
    Rationale: do not mix Modelica build and exec files from different configurations.
    TODO: autoclean (delete) at end
    :param str root_trash_folder: the `.../trash` folder to put trash data in
    :return str: path of the created subfolder"""
    def gen_name(config=None, k=0):
        if config is None:
            return bytes(random.getrandbits(8) for _ in range(3)).hex()
        return format_folder_name(configuration) + f'_{k}'

    get_path = lambda k: path.join(root_trash_folder, gen_name(configuration, k))
    p = get_path(0)
    k = 0
    while path.exists(p):
        k += 1
        p = get_path(k)
    # path `p` does not exist
    os.mkdir(p)
    return p


def check_pipeline_status(output_folder, logger):
    """Check the advancement status of the output_folder.
    :param output_folder: the data folder for a specific configuration
    :return: int, 0 = pipeline has not begun (no folder or empty folder)
                  1 = simulation data is generated
                  2 = pre-processing completed
                 -x = a problem occured at step x
    """
    def check_list(files, checkfun=path.exists):
        paths = map(lambda fp: os.path.join(output_folder, fp), files)
        return list(map(checkfun, paths))

    if not path.exists(output_folder):
        return 0

    step2 = check_list(['X.mat', 'Y.mat'])
    if any(step2) and not all(step2):
        return -2
    elif all(step2):
        return 2

    step1 = check_list(['outputs', 'parameters.txt'])
    if any(step1):
        if not all(step1):
            logger.log('Either `outputs/` or `parameters.txt` is missing')
            return -1
        # Check all output files are present
        output_files = os.listdir(path.join(output_folder, 'outputs'))
        if len(output_files) != N_samples:
            logger.log('`N_samples` do not match with the content of ' + path.basename(output_folder))
            return -1

    return 0


def pipeline(configuration):
    # ------------- PATH & LOGGING CONFIGURATION
    folder_name = format_folder_name(configuration)
    output_folder = path.join(root_output_folder, folder_name)
    # WARNING: don't set this buffer size to 0 (running simulations requires an empty folder)
    logger = Logger(path.join(output_folder, '0-log.log'), print=True, buffer_size=2)

    # ------------- PIPELINE STATUS
    status = check_pipeline_status(output_folder, logger)
    if status < 0:
        logger.log('<ERROR> An error occured while checking ' + folder_name +
                   f', please check by hand [status {status}]. Abortion...')
        logger.flush()
        return
    elif status == 2:
        # Nothing to do...
        return

    # ------------- SIMULATION DATA GENERATION
    if status < 1:
        step1_0D_generation(output_folder, configuration, logger)
    elif status == 1:
        logger.log('Resuming pipeline at step ' + str(status) + ' for ' + folder_name)

    # ------------- PREPROCESSING
    step2_preprocessing(output_folder, logger)

    # ------------- TERMINATE
    logger.log('End of pipeline for ' + folder_name)
    # Flush the buffer
    logger.flush()


def step1_0D_generation(output_folder, configuration, logger):
    def callback(n, step=1000):
        if n % step == 0:
            logger.log(f'Simulation progression: {n}/{N_samples}')

    simlogger = Logger(path.join(output_folder, '1-log_data_generation.log'), buffer_size=200)
    rpm, lvad, artpulse = configuration

    # --- Trash data (modelica build and exec files)
    trash_folder = gen_trash_folder(root_trash_folder, configuration)
    os.chdir(trash_folder)

    # --- Build parameters based on configuration
    params_override = global_0D_params_override.copy()
    params_override.update(build_0D_parameters_from_config(configuration))

    # --- Run simulations
    logger.log('Launching 0D simulations for config ' + path.basename(output_folder))
    uo.runSimulation(
        N=N_samples,
        param_lst=list_parameters,
        output_folder=output_folder,
        file=modelica_file_path,
        LVAD=lvad,
        om_build_settings=om_build_settings,
        override_params=params_override,
        log=simlogger.log,
        om_runtime_args=om_runtime_args,
        callback=callback
    )
    logger.log('Simulations completed for config ' + path.basename(output_folder))

    # --- Save settings & output
    simlogger.flush()
    with open(path.join(output_folder, '1-simulation_settings.json'), 'w') as f:
        json.dump(om_build_settings, f, indent=4)
    with open(path.join(output_folder, '1-params_0D_override.json'), 'w') as f:
        json.dump(global_0D_params_override, f, indent=4)


def step2_preprocessing(output_folder, logger):
    os.chdir(output_folder)
    logger.log('Launching pre-processing for configuration ' + path.basename(output_folder))
    p = subprocess.Popen(['python3', preprocessing_script_path],
                         stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    stdout, stderr = p.communicate()
    with open('2-preprocessing-stdout.log', 'w') as f:
        f.write(stdout.decode() + '\n' + stderr.decode())
    logger.log('Pre-processing has ended')


def gen_training_data(X, Y, logger):
    """
    :param numpy.ndarray X:
    :param numpy.ndarray Y:
    :return:
    """
    os.chdir(root_output_folder)
    if not os.path.exists('dnns'):
        logger.log('Creating folder dnns/')
        os.mkdir('dnns')
    os.chdir('dnns')

    if path.exists('normdata.mat'):
        logger.log('Normalized data already exsits, skipping...')
        return
    else:
        logger.log('Generating normalized data...')
        # --- Data shuffling
        logger.log('Shuffling data...')
        indices = np.random.permutation(X.shape[0])
        X = X[indices]
        Y = Y[indices]

        # --- Data splitting
        logger.log('Splitting data into train and test sets...')
        idx = int(0.95 * X.shape[0])
        Xtrain, Xtest = X[:idx], X[idx:]
        Ytrain, Ytest = Y[:idx], Y[idx:]

        # --- Data normalization
        logger.log('Normalizing data...')
        Xmins = Xtrain.min(axis=0)
        Xmaxs = Xtrain.max(axis=0)
        Xtrain = (Xtrain - Xmins) / (Xmaxs - Xmins)
        Xtest =  (Xtest -  Xmins) / (Xmaxs - Xmins)

        Ymins, Ymaxs = Y.min(0), Y.max(0)
        Ytrain = (Ytrain - Ymins) / (Ymaxs - Ymins)
        Ytest  = (Ytest  - Ymins) / (Ymaxs - Ymins)

        data = {
            'Xmins': Xmins,
            'Xmaxs': Xmaxs,
            'Ymins': Ymins,
            'Ymaxs': Ymaxs,
            'Xtrain': Xtrain,
            'Xtest': Xtest,
            'Ytrain': Ytrain,
            'Ytest': Ytest
        }

        # --- Save data
        logger.log('Saving normalized data...')
        sio.savemat('normdata.mat', mdict=data)


def train_dnn(architecture):
    # Pay attention: 1D arrays stored in .mat files are recovered as 2D array with axis 0 of size 1

    fname = format_dnn_folder_name(architecture)
    if path.exists(fname):
        print(f'Skipping {fname}')
        return
    else:
        os.mkdir(fname)

    data = sio.loadmat('normdata.mat')

    data['Xtrain'] = data['Xtrain']
    data['Ytrain'] = data['Ytrain']

    layers, neurons = architecture
    model = ud.build_keras_model(
        input_shape=(data['Xtrain'].shape[1],),
        noutparams=data['Ytrain'].shape[1],
        n_hlayers=layers, n_neurons=neurons
    )

    start = time.time()
    print(f'Fitting {fname}')
    history = model.fit(data['Xtrain'], data['Ytrain'],
                        epochs=1000, validation_split=0.2, verbose=0)
    print(f'{fname} trained in {time.time() - start:.4} s')

    model.save(path.join(fname, 'model.h5'))
    with open(path.join(fname, 'history.bin'), 'wb') as f:
        pickle.dump(history.history, f)

    plt.plot(history.history['loss'], label='loss')
    plt.plot(history.history['val_loss'], label='val_loss')
    plt.legend()
    plt.savefig(path.join(fname, 'losses.eps'), transparent=False)

    # --- Test
    Ytesthat = model.predict(data['Xtest'])
    Ytest = data['Ytest']
    # Un-normalize
    Ytest = Ytest * (data['Ymaxs'] - data['Ymins']) + data['Ymins']
    Ytesthat = Ytesthat * (data['Ymaxs'] - data['Ymins']) + data['Ymins']
    np.savetxt(path.join(fname, 'Ytest.txt'), Ytest)
    np.savetxt(path.join(fname, 'Ytestpred.txt'), Ytesthat)

    # Simply for compatibility
    normdata = {'parammins': data['Ymins'][0], 'parammaxs': data['Ymaxs'][0]}
    ud.create_and_save_performance_fig(Ytest, Ytesthat, normdata, fname)


# ------------------------------------------------------- #
#                     POST-PROCESSING                     #
# ------------------------------------------------------- #


def concat_logs(root_folder, level='main'):
    """
    :param root_folder:
    :param str level: main = concat main logs (mainlog.log and 0-log.log)
                      all = concat all logs (all files ending in .log)
    :return:
    """
    ls = [path.join(root_folder, 'mainlog.log')]
    filterfun = lambda fname: fname == '0-log.log'
    if level == 'all':
        # TODO: doesn't work
        filterfun = lambda fname: fname.endswith('.log')
        raise NotImplementedError

    for root, folders, files in os.walk(root_folder):
        # Skip trash folder
        if 'trash' in root:
            continue

        logs = [path.join(root, fname) for fname in filter(filterfun, files)]
        ls.extend(logs)

    ls = set(ls)
    df = []
    for fpath in ls:
        print(fpath)
        df.append(pd.read_csv(fpath, header=None, sep=';'))

    df = pd.concat(df)
    df.columns = ['pid', 'time', 'msg']
    df.sort_values('time', inplace=True, ignore_index=True)
    df.time = df.time.apply(datetime.fromtimestamp)
    return df


def save_concat_logs(fpath, df):
    """
    :param str fpath:
    :param pandas.DataFrame df:
    :return:
    """
    df.to_csv(
        fpath,
        sep=';',
        index=False
    )

# ------------------------------------------------------- #
#                        SCRIPTING                        #
# ------------------------------------------------------- #


if __name__ == '__main__':
    main()
