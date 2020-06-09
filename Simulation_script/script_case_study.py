from OMPython import OMCSessionZMQ
import os.path
from os import system


# ========================= INITIALIZATION

output_folder = '/media/maousi/Data/tmp/case_study'

omc = OMCSessionZMQ()


def runcmd(cmd, output=True):
    print(f'Command: {cmd}')
    out = omc.sendExpression(cmd)
    if output:
        print(f'Output: {out}')


runcmd('getVersion()')

model_path = "/media/maousi/Data/Documents/Programmation/git/vad-0d-dnn/modelica/original/Mathcard.mo"
model_name = "Mathcard.Applications.Ursino1998.HMIII.Ursino1998Model_VAD2"

# Load models
runcmd('loadModel(Modelica)')
runcmd(f'loadFile("{model_path}")')
#runcmd(f'instanciateModel("{model_name}")')
runcmd(f'simulate({model_name}, stopTime=30.0, numberOfIntervals=2000,'
       f'simflags=\"-emit_protected\", outputFormat=\"csv\")')


# ========================= DEFINE SIMULATION PARAMETERS

param_artificial_pulse = {
    True: 'LVAD.HMIII_Pulse_Amplitude=2000',
    False: 'LVAD.HMIII_Pulse_Amplitude=0'
}

param_rpm = {
    4000: 'Param_LVAD_RPM=4000',
    5000: 'Param_LVAD_RPM=5000',
    6000: 'Param_LVAD_RPM=6000'
}

param_heart_failure_header = ['Param_LeftVentricle_Emax0', 'Param_LeftVentricle_EmaxRef0',
                              'Param_LeftVentricle_AGain_Emax', 'Param_LeftVentricle_kE']
param_heart_failure = {
    'SHF': [0.2,  0.2,   0.2,   0.011],
    'MHF': [0.8,  0.8,   0.2,   0.013],
    'HH' : [2.95, 2.392, 0.475, 0.014]
}
# Format values: param1_name=value,param2_name=value,...
param_heart_failure = {
    level: ','.join([f'{name}={value}' for name, value in zip(param_heart_failure_header, values)])
    for level, values in param_heart_failure.items()
}

print(param_heart_failure)


# ========================= RUN SIMULATIONS

# Combine all possible configurations of (art_pulse, rpm, hf_level)
# data structure: dictionnary tuple:str
# the key tuple is the human readable configuration
# the value str is the corresponding Modelica parameters
cases = {}
for ap_human, ap_param in param_artificial_pulse.items():
    for rpm_human, rpm_param in param_rpm.items():
        for hf_human, hf_param in param_heart_failure.items():
            key = (ap_human, rpm_human, hf_human)
            params = ','.join([ap_param, rpm_param, hf_param])
            cases[key] = params

print(len(cases))
for k, v in cases.items():
    print(k, v)

output_file_format = 'Ursino1998Model_VAD2_AP_{}_RPM_{}_HF_{}.csv'
command_format = './' + model_name + ' -override={} -r={}'

for human_values, params in cases.items():
    output_file = output_file_format.format(*human_values)
    output_file = os.path.join(output_folder, output_file)
    command = command_format.format(params, output_file)
    print('Running', command)
    system(command)

