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
        * `test_dnn`: all the code related to the DNN testing was moved here
        * `manage_args` to process arguments of [Deep_learning/test_driver.py](Deep_learning/test_driver.py)
    * [Simulation_script/utils_openmodelica.py](Simulation_script/utils_openmodelica.py): 
      new function `runTestSimulation` to generate the DNN test data

## Workflow

Note: some output folder paths have to be configured (e.g. `output_path`),
however, the folder created on the disk might contain the current date at the end.
Those output paths will be used in downstream procedures: be aware that when the 
description refers to the output folder, it implicitly refers to the 
output folder *with* the date as a suffix (e.g. we refer to `output_path` instead of `output_path_MM_dd_YYYY`).

1. Setup simulation by modifying [Simulation_script/setup.py](Simulation_script/setup.py).
    * Make sure that the `output_folder` is in a disk containing enough space (very minimum: 10GB for N=10'000).
      It will be the main folder that contains all the data.
    * Note that `output_folder`'s parent must exist, but but `output_folder` itself !
      This is to ensure that nothing will be overwritten.
    * It is advised to let `DEBUG_MODE = True`, which will make the script `launch_openmodelica.py` *very* verbose 
    so that we can keep track of everything (more details in step 3).
    
1. Create a temporary directory (anywhere, it does not matter and it will contain relatively small files) and `cd` 
   into it. This will contain modelica build files.

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

1. Move `X.mat` and `Y.mat` in the `output_folder` (the one mentionned in steps 3, 4, 5) 
   to group everything together. If you are using Octave, you have to convert the files to binary 
   format first by using `save X.mat -v7` and `save Y.mat -v7`. 

1. `cd` in the `output_folder` and create a directory that will contain data generated
   by the deep neural network, e.g. `mkdir dnn`. At this point, the `output_folder` has the following hierarchy:

		├── dnn/
		├── log.txt
		├── outputs/
		│   ├── Ursino1998Model_VAD2_output_0.mat
		│   ├── Ursino1998Model_VAD2_output_1000.mat
		│   ├── Ursino1998Model_VAD2_output_1001.mat
		│   ├── ... # many other files
		├── parameters.txt
		├── X.mat
		├── Y.mat

1. `cd` into the folder `dnn` and train the DNN: `python <path-to-project>/Deep_learning/test_driver.py train --path ../`.
    
    * It is better to first run `test_driver.py` with the `train` option since no options implies `train` and `test`, 
    which will make the working directory messy (because of the modelica build files). Moreover, testing needs 
    to setup a path, as specified in step 10.
    * The `--path` option is used to find the path to `X.mat` and `Y.mat` files
    * The content of `dnn` will be:
    
			dnn/
			├── coefmaxs.npy
			├── coefmins.npy
			├── DNN_0D_Model.h5
			├── Losses.eps
			├── parammaxs.npy
			├── parammins.npy
			├── Xtest_norm.npy
			└── Ytest_norm.npy

1. Modify `dnn_folder` and `output_folder_DNN_test` in [Simulation_script/setup.py](Simulation_script/setup.py)
   so that they are consistent with the `output_path`. 
    * `dnn_folder` is the folder created at step 8.
    * `output_folder_DNN_test` will contain the DNN test data generated by modelica.
      The folder should not exist yet (but its parent has to exist). Note that the current date will be added at the end of the folder name. 

1. `cd` into a temporary directory (it will contain modelica build files). 
   Test the DNN and simulate with the predicted vs exact parameters with modelica:
   `python <path-to-project>/Deep_learning/test_driver.py test > log_test_dnn.txt &`
    * This will add `DNN_Performance.eps` in the `dnn_folder`
    * This will generate the simulation test data

1. Move `log_test_dnn.txt` to the `output_folder`. The `output_folder` now has the following structure:

		.
		├── dnn  
		│   ├── coefmaxs.npy  
		│   ├── coefmins.npy  
		│   ├── DNN_0D_Model.h5
		│   ├── DNN_Performance.eps
		│   ├── Losses.eps
		│   ├── parammaxs.npy
		│   ├── parammins.npy
		│   ├── Xtest_norm.npy
		│   └── Ytest_norm.npy
		├── dnn_test_2020_04_07/
		│   └── outputs
		│       ├── Ursino1998Model_VAD2_output_0_exact.mat
		│       ├── Ursino1998Model_VAD2_output_0_predicted.mat
		│       ├── ...
		├── log_generation.txt
		├── log_test_dnn.txt
		├── outputs/
		│   ├── Ursino1998Model_VAD2_output_0.mat
		│   ├── Ursino1998Model_VAD2_output_1000.mat
		│   ├── ...
		├── X.mat
		├── Y.mat
