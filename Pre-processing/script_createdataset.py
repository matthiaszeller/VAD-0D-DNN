from buildingspy.io.outputfile import Reader
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os
from natsort import natsorted


def extract_results(varname, reader):
    """reader must be an instance of buildingspy.io.outputfile.Reader"""
    return reader.values(varname)


def time_range(var, t, t_min, t_max):
    indices = np.where(np.logical_and(t >= t_min, t <= t_max))
    return var[indices], t[indices]


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

    isimag = False
    for ak in aks:
        if np.imag(ak) != 0.0:
            isimag = True
            break
    for bk in bks:
        if np.imag(bk) != 0.0:
            isimag = True
            break
    if isimag:
        raise Exception("myyft: some aks, bks are imaginary")

    return coefs, np.real(aks), np.real(bks)


def get_constant_time_steps(t, signal, dt):
    """t in seconds, dt in seconds"""
    df = pd.DataFrame(index = pd.to_datetime(t * 1e9)) # to nanosecond
    df['signal'] = signal
    df = df.resample('{}ms'.format(dt*1000)).mean()
    df.index -= df.index[0]
    newtime = df.index.microseconds * 1e-6 + df.index.seconds
    return newtime.values, df['signal'].values


def perform_fft(signal, t, dt):
    """Equivalent of performmyfftonsignal"""

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

    input_vars = ['SystemicArteries.PC', 'PulmonaryArteries.PC']
    p = '/media/maousi/Data/tmp/simulations_2020_03_21/outputs'
    dircontent = os.listdir(p)
    files = []
    for f in dircontent:
        if f.endswith('.mat'):
            files.append(f)

    # Sort files in natural order, this is CRITICAL since data extracted from files
    # must match the parameters in the file parameters.txt (which obviously uses natural ordering)
    files = natsorted(files)
    X = np.zeros((len(files), len(input_vars), ))

    for f in files:
        reader = Reader(os.path.join(p, f), 'dymola')
        for i in range(len(input_vars)):
            pass
