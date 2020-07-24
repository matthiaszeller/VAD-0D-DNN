"""DO NOT RUN THIS SCRIPT YOURSELF. Use `script_stats_dnn.py`"""

import sys
import os
sys.path.insert(0, os.path.dirname(sys.argv[0]) + '/../Simulation_script/')

import utils_deeplearning as udl
import utils_openmodelica as uo

# Args: n_hlayers, n_neurons, n_aks
args = sys.argv[1:]
n_hlayers, n_neurons, n_aks = [int(e) for e in args]

print(args)

(normdata, (Xtest, Ytest), model) = \
    udl.train_dnn(0.0, files_path='../', save_test_data=False,
    verbose=False, selected_aks=n_aks,
    n_hlayers=n_hlayers, n_neurons=n_neurons)


param1 = uo.Parameter('Param_LeftVentricle_Emax0')
param2 = uo.Parameter('Param_LeftVentricle_EmaxRef0')
param3 = uo.Parameter('Param_LeftVentricle_AGain_Emax')
param4 = uo.Parameter('Param_LeftVentricle_kE')
param_lst = [param1, param2, param3, param4]

udl.test_dnn(model, Xtest, Ytest, normdata, param_lst,
             output_dnn_test=None, modelica_file_path=None,
             dnn_folder='.', runsim=False)


print('Script ended...')
