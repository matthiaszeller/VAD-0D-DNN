"""Run this script directly in the folder containing outputs/ and parameters.txt,
this automatically load the files. Alternatively, you can setup paths in setup_preprocessing.py

WARNING: DO NOT FORGET TO SETUP TIME DISCRETIZATION IN setup_preprocessing.py
"""

from buildingspy.io.outputfile import Reader
import numpy as np
import pandas as pd
import os
from natsort import natsorted
from setup_preprocessing import *
import scipy.io as sio


def extract_results(varname, reader):
    """reader must be an instance of buildingspy.io.outputfile.Reader
    :return: (t, signal)
    """
    return reader.values(varname)


def time_range(var, t, t_min, t_max):
    """
    :return: t, signal
    """
    indices = np.where(np.logical_and(t >= t_min, t <= t_max))
    return t[indices], var[indices]


def myfft(signal):
    n = len(signal) - 1
    h = 2 * np.pi / (n + 1)

    mu = n % 2
    M = int((n - mu) / 2)

    coefs = np.zeros(n + 1, dtype=complex)
    indices = np.arange(0, n + 1)  # 0 to n

    if mu == 1:
        coefs = np.zeros(n + 2, dtype=complex)
        for i in indices:
            #                   don't take indices-1 or i-1 since begin with 0
            val = (np.dot(np.exp(-1j * (indices) * (i - M) * h), signal)) / (n + 1)
            coefs[i + mu] = val

        val = 1 / (2 * (n + 1)) * np.dot(np.power(-1, indices), signal)
        coefs[0] = val
        coefs[-1] = val
    else:
        for i in indices:
            #                   don't take indices-1 or i-1 since begin with 0
            val = (np.dot(np.exp(-1j * (indices) * (i - M) * h), signal)) / (n + 1)
            coefs[i + mu] = val

    if mu == 0:
        midindex = int(np.ceil((n + 1) / 2)) - 1  # -1 since begin with 0
        aks = np.zeros(M + 1, dtype=complex)
        bks = np.zeros(M + 1, dtype=complex)
        # aks[0] = coefs[midindex] * 2
        # bks[0] = 0
        for i in range(0, M + 1):  # 0 to M
            aks[i] = coefs[midindex - i] + coefs[midindex + i]
            bks[i] = 1j * (coefs[midindex + i] - coefs[midindex - i])
            # bks[i+1] = (coefs[midindex + i] - coefs[midindex - i])
    else:
        midindex = int(np.ceil((n + 2) / 2)) - 1  # -1 since begin with 0
        aks = np.zeros(M + 2, dtype=complex)
        bks = np.zeros(M + 2, dtype=complex)
        aks[0] = coefs[midindex] * 2
        bks[0] = 0
        for i in range(0, M + 1):
            aks[i] = coefs[midindex - i] + coefs[midindex + i]
            # bks[i] = 1j * (coefs[midindex+i] - coefs[midindex-i])
            bks[i] = 1j * (coefs[midindex + i] - coefs[midindex - i])
        aks[-1] = 2 * coefs[-1]
        bks[-1] = 0

    # isimag = False
    # for ak in aks:
    #     if np.imag(ak) != 0.0:
    #         isimag = True
    #         break
    # for bk in bks:
    #     if np.imag(bk) != 0.0:
    #         isimag = True
    #         break
    # if isimag:
    #     raise Exception("myyft: some aks, bks are imaginary")

    return coefs, np.real(aks), np.real(bks)


def get_constant_time_steps(t, signal, dt):
    """t in seconds, dt in seconds"""
    df = pd.DataFrame(index = pd.to_datetime(t * 1e9)) # to nanosecond
    df['signal'] = signal
    # Resample to get constant time steps
    # NOTE/ TODO: this is probably not the best way, we should use interpolate only
    # but this creates ValueError: cannot reindex a duplicate axis
    df = df.resample('{}ms'.format(dt*1000)).mean().interpolate('linear')
    # Shift
    df.index -= df.index[0]
    # Create array for time
    newtime = df.index.microseconds * 1e-6 + df.index.seconds
    return newtime.values, df['signal'].values


def perform_fft(signal, t, dt):
    """Equivalent of performmyfftonsignal.m"""

    t, indices = np.unique(t, return_index=True)
    signal = signal[indices]

    t, signal = get_constant_time_steps(t, signal, dt)

    tmin, tmax = t[0], t[-1]
    T = tmax - tmin

    t = t[:-1]
    signal = signal[:-1]
    y, aks, bks = myfft(signal)
    return y, aks, bks, T, t, signal


if __name__ == '__main__':
    # Manage paths
    if path_params is None or path_mats is None:
        path_params = 'parameters.txt'
        path_mats = 'outputs/'

    if not os.path.exists(path_mats) or not os.path.exists(path_params):
        raise Exception('Paths are incorrect.')

    # Create the list of files
    files = []
    for f in os.listdir(path_mats):
        if f.endswith('.mat'):
            files.append(f)

    # Sort files in natural order
    # This is CRITICAL since data extracted from files
    # must match the parameters in the file parameters.txt (which obviously uses natural ordering)
    files = natsorted(files)
    Nfile = len(files)

    # ================ CREATE DNN INPUT DATA X.mat
    for i, f in enumerate(files):
        print(f'Parsing input file {i+1}/{Nfile} {f}')
        reader = Reader(os.path.join(path_mats, f), 'dymola')

        for j, var in enumerate(input_vars):
            # Read the time series corresponding to var in the .mat file
            t, signal = extract_results(var, reader)
            # Get values within tsub_min, tsub_max
            t, signal = time_range(signal, t, tsub_min, tsub_max)
            # Fourier transform
            y, aks, bks, T, t, signal = perform_fft(signal, t, dt)
            # Some bks are always zero
            #print('before reduction, len(bks) =', len(bks))
            #print('index of zero', np.argwhere(bks == 0))
            bks = bks[abs(bks) > 1e-15]
            #print('after reduction, len(bks) =', len(bks))

            # Instanciate null tensor
            if i == 0 and j == 0:
                n = len(aks) + len(bks)
                X = np.zeros((Nfile, len(input_vars), n))

            # Fill in values
            X[i, j, :] = np.append(aks, bks)

    # Save input data
    sio.savemat('X.mat', mdict={'X': X})
    print('Written X.mat')

    # ================= CREATE DNN TARGET DATA Y.mat
    print('Extracting target data')

    Y = pd.read_csv(path_params)
    # Remove column 'n'
    Y.drop('n', axis=1, inplace=True)
    # Extract only the values (no header)
    Y = Y.values
    # Save as .mat file
    sio.savemat('Y.mat', mdict={'Y': Y})
    print('Written Y.mat')

    print('Done!')
