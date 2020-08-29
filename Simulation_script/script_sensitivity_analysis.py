"""Perform 0D model sensitivity analysis. Use -o option to define output folder (default is cwd).

WARNING: you must set the LVAD speed in Mathcard.mo manually (Param_LVAD_RPM)
         before launching this script."""

# ========================================================= #
# ------------------------ IMPORTS ------------------------ #
# ========================================================= #

import utils_openmodelica as uo
from setup import *
import argparse
from os.path import join

# ========================================================= #
# ------------------------- SETUP ------------------------- #
# ========================================================= #

nsamples = 100

# ========================================================= #
# ---------------------- SCRIPT BODY ---------------------- #
# ========================================================= #


def sensitivity_analysis_param(paramlst, id_param_analyzed,
                               nsamples, output_folder):
    """
    :param paramlst: list[uo.Parameter]
    """
    # Reset parameters
    for p in paramlst:
        # Setting parameter value also disables random sampling
        p.setValue( (p.minparam + p.maxparam) / 2 )
    # Activate random sampling for the parameter of interest
    paramlst[id_param_analyzed].randomsampling = True

    # Output folder
    output_folder = join(output_folder, f'P{id_param_analyzed}_{paramlst[id_param_analyzed].name}')

    # Generate simulations
    # The key is to change the default samplingfun argument
    uo.runSimulation(N=nsamples, param_lst=paramlst, output_folder=output_folder,
                     file=file_path, LVAD=SIMULATION_LVAD, samplingfun=uo.sampleParamsEpsilon)


if __name__ == '__main__':
    # Get the output folder from user input
    parser = argparse.ArgumentParser()
    parser.add_argument('-o', '--output_folder', help='default is current working directory',
                        default='.')
    args = parser.parse_args()
    output_folder = args.output_folder

    #output_folder = join(output_folder, datetime.now().strftime('%d-%m-%Y'))

    # Setup parameters for sensitivity analysis
    epsilon = 0.01
    param1 = uo.Parameter("Param_LeftVentricle_Emax0", 0.2, 2.95, epsilon)
    param2 = uo.Parameter("Param_LeftVentricle_EmaxRef0", 0.2, 2.392, epsilon)
    param3 = uo.Parameter("Param_LeftVentricle_AGain_Emax", 0.2, 0.475, epsilon)
    param4 = uo.Parameter("Param_LeftVentricle_kE", 0.011, 0.014, epsilon)
    paramlst = [param1, param2, param3, param4]

    # Perform analysis for each param
    for i in range(4):
        sensitivity_analysis_param(paramlst, i, nsamples, output_folder)
