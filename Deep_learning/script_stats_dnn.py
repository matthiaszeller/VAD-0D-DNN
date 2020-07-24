"""Train different DNN configurations on a given dataset, compare validation losses.

How to run this script:
0) Edit the list `architecture` below
1) Create a folder that will contain DNN statistics for a given dataset
2) Copy X.mat and Y.mat files of the corresponding dataset into that folder
3) cd into that folder, run this script
"""

import tensorflow as tf
from tensorflow import keras
import numpy as np
import sys
sys.path.insert(1, sys.path[0] + '/../Simulation_script/')
import utils_deeplearning as udl
import utils_openmodelica as uo
import os
import random
import shutil
import multiprocessing

# Data structure to specify DNN architectures to train, list of tuples
# (nb hidden layers, nb neurons per hidden layer, number of ak coeffs)
architectures = [
    (1, 8, 40), (1, 16, 40), (1, 32, 40), (1, 64, 40),
    (2, 8, 40), (2, 16, 40), (2, 32, 40), (2, 64, 40),
    (4, 8, 40), (4, 16, 40), (4, 32, 40), (4, 64, 40),
    (6, 8, 40), (6, 16, 40), (6, 32, 40), (6, 64, 40)
]
architectures = [(1, 8, 3)]

routine_script = os.path.dirname(sys.argv[0]) + '/stats_dnn.py'


def format_dir_name(n_hlayers, n_neurons, n_aks):
    return f'{n_hlayers}_hlayers_{n_neurons}_neurons_{n_aks}_aks'


def execute(args):
    n_hlayers, n_neurons, n_aks = args
    print('<DEBUG> Subprocess launched:', n_hlayers, n_neurons, n_aks)
    output_path = format_dir_name(n_hlayers, n_neurons, n_aks)
    # Create the output folder
    os.mkdir(output_path)
    print('<DEBUG> Created folder', output_path)
    # Change working directory to output path
    os.chdir(output_path)
    print('<DEBUG> cd into', output_path)

    # Run external script & train DNN (reminder: this will be launched in subfolder)
    #print('current dir', os.getcwd())
    print('Running rountine...')
    os.system(f'python3 {routine_script} {n_hlayers} {n_neurons} {n_aks} > log.txt')

    # Restore original woking directory - ACTUALLY USELESS
    #os.chdir('..')
    #print('<DEBUG> back to' + os.getcwd())


def main():
    # Manage user errors
    files = os.listdir('.')
    if len(files) != 2 or 'X.mat' not in files or 'Y.mat' not in files:
        raise Exception('You must run the script in a folder containing ONLY 2 files: X.mat and Y.mat')

    # Begin logging
    print('Tensorflow version:', tf.__version__)
    print('\nArchitectures:', architectures)
    print()

    # Multiprocessing: train DNNs in parallel using file `stats_dnn.py`
    print('<DEBUG> Launching multiprocessing...')
    process_pool = multiprocessing.Pool(processes=4)
    process_pool.map(execute, architectures)
    print('<DEBUG> Multiprocessing ended...')


if __name__ == '__main__':
    main()
