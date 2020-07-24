"""DO NOT RUN THIS SCRIPT YOURSELF. Use `script_stats_dnn.py`"""

import sys
import os
sys.path.insert(0, os.path.dirname(sys.argv[0]) + '/../Simulation_script/')

import utils_deeplearning as udl

# Args: n_hlayers, n_neurons, n_aks
args = sys.argv[1:]
n_hlayers, n_neurons, n_aks = [int(e) for e in args]

print(args)

res = udl.train_dnn(0.0, files_path='../', save_test_data=False,
                    verbose=False, selected_aks=n_aks,
                    n_hlayers=n_hlayers, n_neurons=n_neurons)

print('Script ended...')
