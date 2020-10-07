
# TODO: check pump amplitude and other parameters


from buildingspy.io.outputfile import Reader
import os
import matplotlib
from random import sample
import argparse
import pandas as pd
import numpy as np


n_checks, n_successes = 0, 0


def get_values(fpath, varname):
    r = Reader(fpath, 'dymola')
    return r.values(varname)


def pick_files(folder, n=1):
    ls = os.listdir(folder)
    files = sample(ls, n)
    return [
        os.path.join(folder, file)
        for file in files
    ], len(ls)


def get_parameter_value(simulation_file, param_name):
    """
    :param str simulation_file: an output file from an OpenModelica simulation
    :param str param_name: name of the parameter found in `Mathcard.mo`
    :return: float, value of parameter in `simulation_file`
    :raises ValueError:
    """
    t, param = get_values(simulation_file, param_name)
    # Since it is a model parameter, it does not vary and it has only 2 data points
    # at t=0 and t=SimulationTime
    if len(param) != 2:
        raise ValueError('Model parameters should have only two data points')
    if param[0] != param[1]:
        raise ValueError('Unexpected error, initial and final parameter values should be equal.')

    return param[0]


def parse_config_from_folder_name(folder):
    """"""
    if not os.path.exists(folder):
        raise ValueError(f'This folder does not exist: {folder}')

    folder = os.path.basename(folder)
    features = folder.split('_')
    # TODO: parse other config features (bool LVAD, bool AP)
    config = {
        'rpm': int(features[0])
    }
    return config


def get_file_index(file_path):
    """Extract file index from its name.
    Example: Ursino1998Model_VAD2_output_1006.mat -> 1006"""
    return int(file_path.split('_')[-1].split('.')[0])


def check_and_report_error(param_name, value, expected, file, tolerance=None):
    """Check `value` against `expected`. Allow a almost-equality if `tolerance` is not None."""
    global n_checks
    global n_successes
    #print(param_name, value == expected, type(value), type(expected))
    if tolerance is None:
        equal = value == expected
    else:
        equal = abs(value - expected) < tolerance

    n_checks += 1
    if not equal:
        print(f'<ERROR> Inconsistant parameter {param_name} '
              f'(found {value}, expected {expected}, '
              f'difference {abs(value-expected)}) in file {file}')
    else:
        n_successes += 1


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('path', help='root folder containing data generated by `script_pipeline.py`')
    args = parser.parse_args()
    root_output_folder = args.path

    # List folders of datasets
    datasets = [
        os.path.join(root_output_folder, e)
        for e in os.listdir(root_output_folder)
        if 'trash' not in e
    ]
    datasets = filter(os.path.isdir, datasets)
    datasets_size = {}

    nfiles = 3
    model_params = ['Param_LeftVentricle_Emax0', 'Param_LeftVentricle_EmaxRef0',
                    'Param_LeftVentricle_AGain_Emax', 'Param_LeftVentricle_kE']
    tolerance = 1e-12
    print(f'Checking parameter values with tolerance {tolerance}, {nfiles} files per dataset')

    for dataset in datasets:
        dataset_name = os.path.basename(dataset)
        print('Checking dataset', dataset_name, '...')
        outputs = os.path.join(dataset, 'outputs')
        # Parse config from folder name
        config = parse_config_from_folder_name(dataset)
        # Load parameters.txt
        df = pd.read_csv(os.path.join(dataset, 'parameters.txt'), index_col='n')
        # Pick a few output files randomly
        files, n = pick_files(outputs, n=nfiles)
        datasets_size[dataset_name] = n
        for file in files:
            file_index = get_file_index(file)
            # Check parameters
            rpm = get_parameter_value(file, 'Param_LVAD_RPM')
            check_and_report_error('Param_LVAD_RPM', rpm, config["rpm"], file)
            for p in model_params:
                val = get_parameter_value(file, p)
                expected_val = df.loc[file_index, p]
                check_and_report_error(p, val, expected_val, file, tolerance=tolerance)

    print(f'\nChecking done... successes/n_checks = {n_successes}/{n_checks}')
    print()

    a = np.array(list(datasets_size.values()))
    a -= a[0]
    if not np.all(a == 0):
        print(f'<ERROR> All datasets do not have the same size:\n{datasets_size}')


if __name__ == '__main__':
    main()

