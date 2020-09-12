"""This script aims to wrap all code related to:
    - simulation setup
    - simulation data generation
    - data pre-processing
    - DNN training
    - DNN testing
"""

from os import path
import os
import multiprocessing
import random
import sys
sys.path.insert(1, 'Simulation_script/')
from Simulation_script import utils_openmodelica as uo

# TODO: data integrity (or at least procedure completeness) check


# ========================================================= #
# ------------------------- SETUP ------------------------- #
# ========================================================= #

# ---------- GENERAL SETUP
# Root folder containing all data
root_output_folder = '/media/maousi/Data/tmp/pipelining'
root_trash_folder = path.join(root_output_folder, 'trash')
# Choose this value to get ~ 100% of CPU usage during simulation data generation
N_processes = 3

# All pump configurations (pump speed, LVAD, artificial Pulse)
configurations = [
    (rpm, True, True)
    for rpm in [4000]
]

# ---------- SIMULATION SETUP
# The number of simulations in a single dataset
N_samples = 10
# Whether to use the 0D model with/without LVAD
LVAD_simulation = True
# Simulation time, number of intervals, etc...
om_simulation_settings = {
    'stopTime': 30.0,
    'numberOfIntervals': 2000,
    'outputFormat': '"mat"',
    'simflags': '"-emit_protected"',
}
# 0D Model parameters overriding
global_0D_params_override = {

}
# The path to `Mathcard.mo`
modelica_file_path = '/media/maousi/Data/Documents/Programmation/git/vad-0d-dnn/modelica/original/Mathcard.mo'

# --- Do not edit
param1 = uo.Parameter("Param_LeftVentricle_Emax0", 0.2, 2.95)
param2 = uo.Parameter("Param_LeftVentricle_EmaxRef0", 0.2, 2.392)
param3 = uo.Parameter("Param_LeftVentricle_AGain_Emax", 0.2, 0.475)
param4 = uo.Parameter("Param_LeftVentricle_kE", 0.011, 0.014)
list_parameters = [param1, param2, param3, param4]
# ---

# ---------- SETUP CHECKUP
if path.isfile(root_output_folder):
    raise ValueError('Invalid root_output_folder')
elif not path.exists(root_output_folder):
    os.mkdir(root_output_folder)
else:
    # TODO: determine behaviour in case folder exists,
    # one should check data integrity and procedure completeness
    pass


# ========================================================= #
# ---------------------- SCRIPT BODY ---------------------- #
# ========================================================= #

def format_folder_name(configuration):
    rpm, lvad, artpulse = configuration
    lvad = '_LVAD' if lvad else ''
    artpulse = '_AP' if artpulse else ''
    return f'{rpm}{lvad}{artpulse}'


def build_0D_parameters_from_config(configuration):
    rpm, lvad, artpulse = configuration
    dic = {}
    if not artpulse:
        dic['HMIII_Pulse_Amplitude'] = 0
    return dic


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
                   'LVAD_simulation', 'om_simulation_settings',
                   'global_0D_params_override', 'configurations']

    for var in prompt_vars:
        display(var, globals())

    i = input('\nDo you want to continue ? [y/N] ')
    if i != 'y':
        exit()


def main():
    prompt_initial_validate()

    # Prepare
    if not path.exists(root_trash_folder):
        os.mkdir(root_trash_folder)

    # Go!
    try:
        ps_pool = multiprocessing.Pool(processes=N_processes)
        ps_pool.map(pipeline, configurations)
    except Exception as e:
        print(e)
        print('exiting after error...')
        exit()


def gen_trash_folder(root_trash_folder, configuration=None):
    """Create a folder with a unique name to perform temporary computations.
    Rationale: do not mix Modelica build and exec files from different configurations.
    TODO: autoclean (delete) at end
    :param str root_trash_folder: the `.../trash` folder to put trash data in
    :return str: path of the created folder"""
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


def pipeline(configuration):
    # ------------- SIMULATION DATA GENERATION
    output_folder = path.join(root_output_folder, format_folder_name(configuration))
    # --- Trash data (modelica build and exec files)
    trash_folder = gen_trash_folder(root_trash_folder, configuration)
    os.chdir(trash_folder)
    # --- Build parameters based on configuration
    params_override = global_0D_params_override.copy()
    params_override.update(build_0D_parameters_from_config(configuration))

    # --- Run simulations
    uo.runSimulation(
        N=N_samples,
        param_lst=list_parameters,
        output_folder=output_folder,
        file=modelica_file_path,
        LVAD=LVAD_simulation,
        om_sim_settings=om_simulation_settings,
        override_params=params_override
    )



if __name__ == '__main__':
    main()

