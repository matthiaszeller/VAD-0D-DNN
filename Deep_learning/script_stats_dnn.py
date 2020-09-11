"""Train different DNN configurations on a given dataset, compare losses.

=== Overview:
    * Train several DNNs in parallel to reduce computation time
    * Separately specify each number of layers, neurons, coefficients and combine them all automatically
    * Don't worry about running this script several times in the same folder,
        - nothing gets deleted
        - no recomputation of pre-existing data (detected if corresponding folder already exists),
          those are simply ignored
    * 2 modes: no arguments = train DNNs, --parse = analyze results

=== HOW TO? COMPUTE DIFFERENT ARCHITECTURES / DNN INPUT SIZES
0) Edit the dict `combinations` and the variable `n_processes` below
    -> Note: the number of simultaneously-launched processes (`n_processes`) should be
             chosen in order to have 100% CPU usage (not necessarily the actual number of
             cores since tensorflow distributes computations evenly among all cores)
1) Create a folder that will contain DNN statistics for a given dataset
2) Copy X.mat and Y.mat files of the corresponding dataset into that folder
3) cd into that folder, run this script:
    $ python3 script_stats_dnn.py
4) Eventually, you should have a similar root folder contents:
    .
    ├── 1_hlayers_16_neurons_100_aks/
    ├── 1_hlayers_16_neurons_10_aks/
    ├── [...]
    ├── 6_hlayers_8_neurons_40_aks/
    ├── 6_hlayers_8_neurons_60_aks/
    ├── X.mat
    └── Y.mat

=== HOW TO? EXTRACT, AGGREGATE & SAVE RESULTS
1) Re-run this script with the --analyze option: this creates a .csv and figures
   in the (created) subfolder results.
"""

# ======================================================= #
# ----------------------- IMPORTS ----------------------- #
# ======================================================= #

import sys
sys.path.insert(1, sys.path[0] + '/../Simulation_script/')
import os
import multiprocessing
import argparse
import pickle
import pandas as pd
from time import time

# Disable tensorflow debugging output, only show warnings/errors
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '1'
import tensorflow as tf

# ========================================================= #
# ------------------------- SETUP ------------------------- #
# ========================================================= #

# Choose number of simultaneously-launched processes
n_processes = 6

# [Do not edit] Path of externally-launched script for multiprocessing
routine_script = os.path.dirname(sys.argv[0]) + '/stats_dnn.py'

# Data structure to specify DNN architectures to train, list of tuples
# (nb hidden layers, nb neurons per hidden layer, number of ak coeffs)

# architectures = [
#     (1, 8, 40), (1, 16, 40), (1, 32, 40), (1, 64, 40),
#     (2, 8, 40), (2, 16, 40), (2, 32, 40), (2, 64, 40),
#     (4, 8, 40), (4, 16, 40), (4, 32, 40), (4, 64, 40),
#     (6, 8, 40), (6, 16, 40), (6, 32, 40), (6, 64, 40)
# ]
# architectures = [(1, 8, 40), (4, 16, 40)]

# Define number of each variable (layers, neurons, coeffs)
combinations = {
    'hlayers': [1],#[1, 2, 3, 4, 5, 6],
    'neurons': [8, 16, 32, 64, 128],
    'aks':     [50]
               #[10, 20, 40, 60, 100]
}

# ========================================================= #
# ---------------------- SCRIPT BODY ---------------------- #
# ========================================================= #


def architectures():
    """Architectures generator function: combine every item of each variable (layers,neurons,coeffs)"""
    for aks in combinations['aks']:
        for neurons in combinations['neurons']:
            for hlayers in combinations['hlayers']:
                yield (hlayers, neurons, aks)


def format_dir_name(n_hlayers, n_neurons, n_aks):
    return f'{n_hlayers}_hlayers_{n_neurons}_neurons_{n_aks}_aks'


def execute(args):
    """Called by processes in the pool of processes managed by multiprocessing package."""
    n_hlayers, n_neurons, n_aks = args
    msg = f'[{n_hlayers}, {n_neurons}, {n_aks}]'

    print('<Multiprocessing>', msg, 'Subprocess launched')
    output_path = format_dir_name(n_hlayers, n_neurons, n_aks)
    # Create the output folder
    if os.path.exists(output_path):
        # In case folder exists, computations have already been done
        # => abort the current task and do next one
        print('\n', '='*20, '\n<ERROR> <Multiprocessing>', msg,
              'The folder ' + output_path + ' already exists ! Skipping...\n',
              '='*20, '\n')
        return

    os.mkdir(output_path)
    print('<Multiprocessing>', msg, 'Created folder', output_path)
    # Change working directory to output path
    os.chdir(output_path)
    print('<Multiprocessing>', msg, 'cd into', output_path)

    # Run external script & train DNN (reminder: this will be launched in subfolder)
    #print('current dir', os.getcwd())
    print('<Multiprocessing>', msg, 'Running external script from', os.getcwd())
    start = time()
    os.system(f'python3 {routine_script} {n_hlayers} {n_neurons} {n_aks} > log.txt')
    print('Elapsed time:', time() - start)

    # Restore original woking directory -> required for next task to work
    os.chdir('..')
    print('<Multiprocessing>', msg, 'back to' + os.getcwd())


def manage_tensorflow_gpu(gpu_forced):
    print('Tensorflow version:', tf.__version__)
    print()

    # Display recognized GPUs
    gpus = tf.config.experimental.list_physical_devices('GPU')
    N = len(gpus)
    if N == 0:
        if gpu_forced is True:
            print('<ERROR> No GPU available, but launched with --gpu. Abort.')
            exit()
        print('<WARNING> No GPU available, computations will take place on CPU')
    else:
        print(f'{N} GPU{"s" if N > 1 else ""} recognized:\n{gpus}')


def main():
    # === MODE: ANALYZE RESULTS
    # User arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('--analyze', action='store_true',
                        help='extract results and save results')
    # Note that if computations will automatically take place on the GPU if possible,
    # but the user might forget to activate to correct environment
    parser.add_argument('--gpu', action='store_true',
                        help='force computations on GPU, throw error if not available')
    args = parser.parse_args()
    if args.analyze:
        analyze_results()
        return

    # === MODE: TRAIN DNNS
    # Manage user errors
    files = os.listdir('.')
    if 'X.mat' not in files or 'Y.mat' not in files:
        raise Exception('You must run the script in a folder containing X.mat and Y.mat files.')

    manage_tensorflow_gpu(args.gpu)

    # Begin logging
    print('\nArchitectures:', architectures)
    print()

    # Multiprocessing: train DNNs in parallel using file `stats_dnn.py`
    print('<DEBUG> Launching multiprocessing...')
    process_pool = multiprocessing.Pool(processes=n_processes)
    process_pool.map(execute, architectures())
    print('<DEBUG> Multiprocessing ended...')


def analyze_results():
    """Read and save results stored in the current working directory (cwd)."""
    # === ROUTINES
    def get_list(rootdir, fileformat):
        """
        :param fileformat: example: `.txt`
        """
        ls = []
        for dirname, subdirlist, filelist in os.walk(rootdir):
            for file in filelist:
                if file.endswith(fileformat):
                    ls.append(os.path.join(dirname, file))

        return ls

    def parse_config(filepath):
        parent_dir = os.path.basename(os.path.dirname(filepath))
        ls = parent_dir.split('_')
        res = {
            'hlayers': 0,
            'neurons': 2,
            'aks': 4
        }
        return {key: int(ls[index]) for key, index in res.items()}

    def load_losses(filepath):
        with open(filepath, 'rb') as f:
            data = pickle.load(f)

        # We extract the last item of every list
        return {key: lst[-1] for key, lst in data.items()}

    def process_stats_data(rootdir):
        files = get_list(rootdir, 'history.bin')
        data = []
        N = len(files)
        for i,f in enumerate(files):
            info = parse_config(f)
            print(f'Parsing {i+1}/{N}, {info}')
            info.update(load_losses(f))
            data.append(info)
        return pd.DataFrame(data).sort_values(['hlayers', 'neurons', 'aks'])

    # === FUNCTION BODY
    df = process_stats_data(rootdir='.')
    df.to_csv('results.csv', index=False)
    print('\nOutput written to results.csv ...\n')


if __name__ == '__main__':
    main()
