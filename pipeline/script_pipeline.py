"""This script wraps all the code related to:
    - simulation setup
    - simulation data generation
    - data pre-processing
    - DNN training
    - DNN testing - TODO: implement test data generation

When running the script, you are first prompted to continue.
Advice:
    - to do some tests, set `N_samples` to a very low value (but this will trigger error for DNN training)
    - read the configuration displayed when you run the script
    - use `script_stats_dnn.py` to find the best architectures before running this script
    - don't run more than `N_processes` configurations at once, otherwise somehow the HDD gets busy ?
        # TODO: investigate the problem with HDD
    - you can use `script_check_sim_data.py` to check consistency of simulation data

General script notions:
    - A dataset is made of `N_samples` simulation files, preprocessing files and DNN data
    - A dataset is characterized by its global simulation settings: pump speed (RPM),
      whether the LVAD is used, and whether artificial pulse (AP) is activated
    - The word "pipeline" used in this script refers to generating data for a single dataset
    - This script runs `N_processes` and assigns 1 pipeline to each process
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
    4. repeat from 1. and set the other configurations you want (don't change N_samples)
    5. once all simulation data is computed, you can now train additional DNN configurations,
       set the `configurations` variable accordingly
    5. re-run the script - it trains the additional DNNs (0D data is not computed again)
"""


# ======================================================= #
# ----------------------- IMPORTS ----------------------- #
# ======================================================= #

from os import path
import os
import multiprocessing
import random
import json
import subprocess
import sys
import pandas as pd
import argparse
from datetime import datetime

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
    (rpm, True, True, [(2, 64, 50)])
    for rpm in range(4600, 6100, 100)
    #[4000, 4100, 4200, 4300, 4400, 4500]
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


# ========================================================= #
# ---------------------- SCRIPT BODY ---------------------- #
# ========================================================= #


# --------------------------------------------------------- #
#                     DATASET SETTINGS                      #
# --------------------------------------------------------- #


def format_folder_name(configuration):
    rpm, lvad, artpulse, _ = configuration
    lvad = '_LVAD' if lvad else ''
    artpulse = '_AP' if artpulse else ''
    return f'{rpm}{lvad}{artpulse}'


def build_0D_parameters_from_config(configuration):
    rpm, lvad, artpulse, _ = configuration
    dic = {}
    # --- Pump speed
    dic['Param_LVAD_RPM'] = rpm
    # --- Artificial Pulse
    if not artpulse:
        dic['HMIII_Pulse_Amplitude'] = 0

    return dic


def format_dnn_folder_name(architecture):
    layers, neurons, aks = architecture
    return f'dnn_{layers}_layers_{neurons}_neurons_{aks}_aks'


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
    print('\tRPM,   LVAD, AP,   Architectures')
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
                        help='concatenate logs, write in file and exit',
                        type=str, default=None)
    parser.add_argument('-p', '--path',
                        help='override the root folder, mainly used for debugging',
                        type=str, default=None)
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
    prompt_initial_validate()

    # ------------- PREPARATION
    mainlogger = Logger(path.join(root_output_folder, 'mainlog.log'),
                        print=True, buffer_size=0)

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

    # ------------- TERMINATE
    server_handler.kill_sessions()

    # Concatenate log files
    mainlogger.log('Concatenating log files...')
    df = concat_logs(root_output_folder)
    save_concat_logs(path.join(root_output_folder, 'logs.log'), df)

    # Remove trash files
    if path.exists(root_trash_folder):
        mainlogger.log('Deleting trash data...')
        os.rmdir(root_trash_folder)

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
    :return: int, 0 = output_folder does not exist
                  1 = simulation data is generated
                  2 = pre-processing completed
                 -x = a problem occured at step x
    """
    def check_list(files, checkfun=path.exists):
        fdic = {
            f: path.join(output_folder, f) for f in files
        }
        return {
            f: checkfun(fpath) for f, fpath in fdic.items()
        }

    if not path.exists(output_folder):
        return 0

    step = 1
    step1 = check_list(['outputs', 'parameters.txt'])
    if any(step1) and not all(step1):
        logger.log('One file is missing')
        return -step
    # Check all output files are present
    output_files = os.listdir(path.join(output_folder, 'outputs'))
    if len(output_files) != N_samples:
        logger.log('`N_samples` do not match with the content of ' + path.basename(output_folder))
        return -step

    step = 2
    step2 = check_list(['X.mat', 'Y.mat'])
    if any(step2) and not all(step2):
        logger.log('One pre-processing file is missing')
        return -step

    return step


def pipeline(configuration):
    # ------------- PATH & LOGGING CONFIGURATION
    folder_name = format_folder_name(configuration)
    output_folder = path.join(root_output_folder, folder_name)
    logger = Logger(path.join(output_folder, '0-log.log'), print=True, buffer_size=0)

    # ------------- PIPELINE STATUS
    status = check_pipeline_status(output_folder, logger)
    if status < 0:
        logger.log('<ERROR> An error occured while checking ' + folder_name +
                   f', please check by hand [status {status}]')
        logger.flush()
        return

    # ------------- SIMULATION DATA GENERATION
    if status < 1:
        step1_0D_generation(output_folder, configuration, logger)
    elif status == 1:
        logger.log('Resuming pipeline at step ' + str(status) + ' for ' + folder_name)

    # ------------- PREPROCESSING
    if status < 2:
        step2_preprocessing(output_folder, logger)
    elif status == 2:
        logger.log('Resuming pipeline at step ' + str(status) + ' for ' + folder_name)

    # ------------- DNN TRAINING
    rpm, lvad, artpulse, architectures = configuration
    step3_train_dnns(output_folder, architectures, logger)

    # ------------- TERMINATE
    # Flush the buffer
    logger.flush()


def step1_0D_generation(output_folder, configuration, logger):
    def callback(n, step=1000):
        if n % step == 0:
            logger.log(f'Simulation progression: {n}/{N_samples}')

    simlogger = Logger(path.join(output_folder, '1-log_data_generation.log'), buffer_size=200)
    rpm, lvad, artpulse, architectures = configuration

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


def step3_train_dnns(output_folder, architectures, logger):
    logger.log('Launching DNN training for configuration ' + path.basename(output_folder))
    n = len(architectures)

    for i, (hlayers, neurons, aks) in enumerate(architectures):
        folder_name = format_dnn_folder_name((hlayers, neurons, aks))
        dnn_data_folder = path.join(output_folder, folder_name)
        if path.exists(dnn_data_folder):
            logger.log(f'Architecture {folder_name} already exists, skipping...')
            continue
        logger.log(f'Training architecture {i+1}/{n} - {folder_name}')
        os.mkdir(dnn_data_folder)
        os.chdir(dnn_data_folder)

        # Redirect stdout to a log file for the training/testing procedure
        dnnlogger = Logger(path.join(dnn_data_folder, 'log_dnn.log'))
        old_stdout = sys.stdout
        sys.stdout = dnnlogger

        (normdata, (Xtest, Ytest), model) = ud.train_dnn(
            perccoef=0,
            selected_aks=aks,
            n_hlayers=hlayers,
            n_neurons=neurons,
            files_path='../'
        )

        ud.test_dnn(
            model=model,
            Xtest=Xtest,
            Ytest=Ytest,
            normdata=normdata,
            param_lst=list_parameters,
            runsim=False,
            dnn_folder='.',
            output_dnn_test=None,
            modelica_file_path=None,
        )

        dnnlogger.flush()
        sys.stdout = old_stdout

    logger.log('DNN training completed...')


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
