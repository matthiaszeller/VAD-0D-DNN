# VAD-0D-DNN

## About

Semester project with Jean Bonnemain and Simone Deparis. 
All the code was provided by [Jean Bonnemain](https://bitbucket.org/%7Bbf4fa678-951a-491f-a0d2-7efb1e5d3c9c%7D/). 
The project aims to simulate the cardiovascular system - optionally with a left ventricular assist device (LVAD) - 
and train a deep neural network to estimate parameters reflecting heart disfunctions.
My contribution is to reimplement the code for a new cardiac device, and train a deep neural network 
to estimate parameter values that reflect left ventricular recovery. 

## Overview 

1. Implementation of the 0D model for the cardiovascular system, with LVAD. Performed with Modelica. 
1. Dataset generation by the 0D model 
1. Preprocessing
1. Neural network training and testing
1. Postprocessing

## Brief changelog

1. **HQ curve**: Update the modelica code for Heart Mate III since it is currently implemented for Heart Mate II.
	
    * Fit quadratic polynomials to the HQ curves in the HMIII manual
    * Duplicate the `VAD` model to `VAD2` and modify the expression of `Q` as a function of `dP`
    * Duplicate the main model that instanciates the `VAD` model and instanciate `VAD2` instead
    * Check the consistency of curves for medium heart failure (MHF) / severe heart failure (SHF) 
        and in context of baroreceptor regulation / no regulation. 
        Absence of regulation is simulated by setting the following parameters to zero: 
        `Param_LeftVentricle_Emax0`, `Param_LeftVentricle_EmaxRef0`, 
        `Param_LeftVentricle_AGain_Emax`, `Param_LeftVentricle_kE`

2. **Simulation data generation**: Use the model with `VAD2` to generate the simulation data. 

    * [Simulation_script/utils_openmodelica.py](Simulation_script/utils_openmodelica.py):
    change the data generation procedure, new `runSimulation` function replacing
    `launchSimulation`. Idea: we first build the modelica file to create an executable,
    then we run the executable with `-override`. This significantly improves speed of data generation.
    * [Simulation_script/setup.py](Simulation_script/setup.py): put all parameters here, 
    as the number of samples, the paths, ...
    * Folder/files generated at end of simulation (the only difference is that the content of `models/` is 
    replaced by `parameters.txt`): 
        * `output_folder` is created
        * `output_folder/outputs` is created
        * Output files: As many as `numberofsamples` .mat files are generated in `output_folder/outputs/`
        * Input file: a single `output_folder/parameters.txt` file in .csv format

3. **Preprocessing**: create a dataset (`.mat` files) for the neural network. The inputs are the 
   PAS and PAP curves (Fourier coefficients) and the outputs are the parameters
   reflecting cardiac dysfunction.
    * New [Pre-processing/setup.m](Pre-processing/setup.m) file regrouping some parameters
    (`tsub_min`, `tsub_max`, `dt`, paths...)
    * Adaptation of [Pre-processing/test_createdataset.m](Pre-processing/test_createdataset.m):
        * Convert output data from .csv to .mat format (`parameters.txt` -> `Y.mat`)
    * New functions:
        * [Pre-processing/contains.m](Pre-processing/contains.m): this function is 
        *useless in Matlab* but guarantees compatibility with Octave
        * [Pre-processing/sort_nat.m](Pre-processing/sort_nat.m):
        to sort an Matlab array of files in the natural order. For instance, tansform 
        `{'f1.txt'  'f10.txt'  'f2.txt'}` to `{'f1.txt'  'f2.txt'  'f10.txt'}`

4. **Deep learning**: load the `.mat` files, train and test the neural network.
    * [Deep_learning/test_driver.py](Deep_learning/test_driver.py): accepts arguments (`help`, `train`, `test`), 
    main file content moved in [Deep_learning/utils_deeplearning.py](Deep_learning/utils_deeplearning.py).
    * [Deep_learning/utils_deeplearning.py](Deep_learning/utils_deeplearning.py): new functions:
        * `train_dnn`: all the code related to the DNN training was moved here
        * `manage_args` to process arguments of [Deep_learning/test_driver.py](Deep_learning/test_driver.py)

## Workflow

1. Setup simulation by modifying [Simulation_script/setup.py](Simulation_script/setup.py).
   Make sure that the `output_folder` is in a disk containing enough space (very minimum: 10GB for N=10'000)

1. Create a temporary directory (anywhere, it does not matter) and `cd` into it . 
   This will contain modelica build files.

1. Run `python <path-to-project>/Simulation_script/launch_openmodelica.py > log.txt &` and replace
   `<path-to-project>` accordingly. 
    * It is advised to redirect the stdoutput with `> log.txt` since the script is *very verbose*
    * Build files are created in the (temporary) working directory while
   output files are created in `output_folder`.

1. Once the simulation is finished, you can move `log.txt` to `output_folder`. If one needs to check something
   later, one can find it in the log file.

1. Modify [Pre-processing/setup.m](Pre-processing/setup.m) to setup the paths. 
    * `pathmats` should target the subdirectory `outputs/` located in the directory created
      by [Simulation_script/launch_openmodelica.py](Simulation_script/launch_openmodelica.py)
    * `pathparams` should target the file `parameters.txt` located in the directory created
      by [Simulation_script/launch_openmodelica.py](Simulation_script/launch_openmodelica.py)

1. Open Matlab/Octave and run [Pre-processing/test_createdataset.m](Pre-processing/test_createdataset.m).
   This will generate `X.mat` and `Y.mat` files  (in Matlab working directory).

1. (Optional) Move `X.mat` and `Y.mat` in the `output_folder` (the one mentionned in steps 3, 4, 5) 
   to group everything together.

1. (Optional) `cd` in the `output_folder` and create a directory that will contain data generated
   by the deep neural network, e.g. `mkdir dnn`. `cd` into this directory.
   
1. Run `python <path-to-project>/Deep_learning/test_driver.py`. If `X.mat` and `Y.mat` files are not 
   in the working directory, use the `--path` argument. For instance, if you followed steps 7 and 8,
   run `python <path-to-project>/Deep_learning/test_driver.py --path ../` to 
   load the matrix files from the parent directory.  

## Data

#### Modelica

* [Mathcard.mo](modelica/Mathcard.mo): modelica file containing the suitable libraries to model the cardiovascular system

#### Notebooks

1. [HQ curve HMIII](notebooks/N1-HQ-curve-HMIII.ipynb): fit quadratic polynomials of the pressure head-flow curves for different RPMs of the Heart Mate III
