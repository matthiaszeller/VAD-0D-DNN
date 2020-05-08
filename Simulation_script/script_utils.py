import argparse
from buildingspy.io.outputfile import Reader
from os.path import join
import pandas as pd
import matplotlib.pyplot as plt

# Useful variables for Interpreter debugging
p = '/media/maousi/Data/tmp/simulation_LVAD_RPM6000_Pulse_T30_N2000_2020_04_26/outputs'
f = 'Ursino1998Model_VAD2_output_1000.mat'
fp = join(p, f)


def extract(file, var):
    reader = Reader(file, 'dymola')
    t, signal = reader.values(var)
    return t, signal


def interact(file, variable, output = 'df'):
    if file is not None and variable is not None:
        t, signal = extract(file, variable)

        if output == 'df':
            df = pd.DataFrame(index=t)
            df[variable] = signal
            return df

        data = f't,{variable}\n'
        for i in range(len(t)):
            data += f'{t[i]},{signal[i]}\n'

        if output is None:
            print(data, end='')
        else:
            with open(output, 'w') as f:
                f.write(data)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--file', help='Path to .mat file')
    parser.add_argument('-v', '--variable', help='Name of the variable to extract')
    parser.add_argument('-o', '--output', help='Path of output file (if not specified, print to stdout)')
    args = parser.parse_args()

    interact(args.file, args.variable, args.output)

if __name__ == '__main__':
    main()
