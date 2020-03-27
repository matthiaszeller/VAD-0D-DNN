
# ======================= ABOUT
# Launch this file to run the simulations.
# Note that many files will be created in the working directory,
# and the output is very verbose (in case DEBUG_MODE is True).
# Advice: run this script and re-direct the output using `>` to keep track of everything
# Example: run `python3 launch_openmodelica.py > log.txt`
# Full procedure:
#   1. In `setup.py`, modify:
#       1.1. `file_path` to match the path of the modelica file
#       1.2. `output_folder` as desired, outputs in .mat format will be written there
#            Note that the folder should not exist and will be created
#   2. Create a temporary directory to store build files (path does not matter)
#   3. `cd` into this temporary directory
#   3. Run `python3 <path-to-project>/Simulation_script/launch_openmodelica.py > log.txt`
#   4. Keep track of data generation with `tail log.txt`
#      WARNING: check that simulations actually run, since a mismatch with the paths
#               may happen (or the output folder may already exist)

# ======================= IMPORTS

import utils_openmodelica as uo
import datetime
# File paths
from setup import *

# ======================= IMPLEMENTATION

today = datetime.datetime.now()
output_folder += "_" + today.strftime("%Y") + '_' + today.strftime("%m") + '_' + today.strftime("%d")

q("<INFO> The output folder is "+output_folder)

# list of parameters which we vary
param1 = uo.Parameter("Param_LeftVentricle_Emax0", 0.2, 2.95)
param2 = uo.Parameter("Param_LeftVentricle_EmaxRef0", 0.2, 2.392)
param3 = uo.Parameter("Param_LeftVentricle_AGain_Emax", 0.2, 0.475)
param4 = uo.Parameter("Param_LeftVentricle_kE", 0.011, 0.014)

listParameters = [param1, param2, param3, param4]

uo.runSimulation(numberofsamples, listParameters, output_folder, file_path, True)
