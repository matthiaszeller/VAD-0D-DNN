import argparse
from buildingspy.io.outputfile import Reader


def extract(args):
    reader = Reader(args.file, 'dymola')
    t, signal = reader.values(args.variable)
    return t, signal


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--file', help='Path to .mat file')
    parser.add_argument('-v', '--variable', help='Name of the variable to extract')
    parser.add_argument('-o', '--output', help='Path of output file (if not specified, print to stdout)')
    args = parser.parse_args()

    if args.file is not None and args.variable is not None:
        t, signal = extract(args)

        data = f't,{args.variable}\n'
        for i in range(len(t)):
            data += f'{t[i]},{signal[i]}\n'

        if args.output is None:
            print(data, end='')
        else:
            with open(args.output, 'w') as f:
                f.write(data)



if __name__ == '__main__':
    main()
