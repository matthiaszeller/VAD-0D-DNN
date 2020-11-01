# VAD-0D-DNN

**TODO: description**

# Getting started

## Prerequisites



# Project Overview

The project is subdivided in the following sections:

* 0D Model Implementation: `modelica`, 0D Model of the cardiovascular system implemented in the OpenModelica language

* 0D Model Simulation: `Simulation_script`, Python code gluing with the model implemented with OpenModelica language 
    * Generate a dataset by solving the model with different parameter values
    * Perform a sensitivity analysis: fix 3 heart-failure parameters and vary the 4th one
    * Case study: simulate different scenarios by varying heart failure severity, pump speed

* Preprocessing: `Pre-processing`, Python and Matlab code 
    * Extract systemic arterial pressure and pulmonary arterial pressure from simulation data
    * Compute Fourier coefficients
    * Generate DNN input and output files

* DNN Training and Testing: `Deep_learning`, Python (Keras Library)
    * DNN selection: train different DNN architectures with different input sizes to select the best DNN
    * Splits dataset into training / test set 
    
* Post-processing: `Post-processing_Code`, Matlab
    * Compute hemodynamic quantities on the test set

* Notebooks: `notebooks`, Jupyter Notebooks (Python)
    * Data analysis and visualization
    * Generate figures and tables for the article
    * Results aggregation
    * Almost each step of the pipeline has one or more dedicated notebook
    
* Project pipelining: `pipelining` Python code gluing data generation, preprocessing, DNN training
    * Automatic way to generate datasets with different pump configurations, merge them, train DNNs


Note that the contents of `pipelining` was used to generate data for the article.


# Workflow

The contents of `Simulation_script` can generate a dataset under a single pump configuration, according to the 
content of the `Mathcard.mo` file. We first each step of the project pipeline and how to manually run it. 
We then describe the automatic method, which wraps all the manual steps (up to DNN training) in a single Python script.

## Manual pipeline

![dataflow](res/diagram/data_flow.bmp)

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

1. **NOTE: you can use Python for the pre-processing step. Go to step 9 to use Python**. Note that the Python script is slower, but it is useful in
some cases when the Matlab code doesn't work. Matlab works for simulation time = 20s, dt = 0.04s. Reported issue: simulation time 30s, dt = 0.015.

1. **Preprocessing with MATLAB**. Modify [Pre-processing/setup.m](Pre-processing/setup.m) to setup the paths. 
    * `pathmats` should target the subdirectory `outputs/` located in the directory created
      by [Simulation_script/launch_openmodelica.py](Simulation_script/launch_openmodelica.py)
    * `pathparams` should target the file `parameters.txt` located in the directory created
      by [Simulation_script/launch_openmodelica.py](Simulation_script/launch_openmodelica.py)

1. Open Matlab and run [Pre-processing/test_createdataset.m](Pre-processing/test_createdataset.m).
   This will generate `X.mat` and `Y.mat` files  (in Matlab working directory).

1. Move `X.mat` and `Y.mat` in the `output_folder` (the one mentionned in steps 3, 4, 5) 
   to group everything together. If you are using Octave, you have to convert the files to binary 
   format first by using `save X.mat -v7` and `save Y.mat -v7`. 

1. **Preprocessing with Python** (skip this step if you used Matlab). 
    1. Setup time discretization in [Pre-processing/setup_preprocessing.py](Pre-processing/setup_preprocessing.py)
    1. `cd` in the `output_folder`
    1. Run [Pre-processing/script_createdataset.py](Pre-processing/script_createdataset.py).
   It will automatically find `parameters.txt` and `outputs/`
    1. It creates `X.mat` and `Y.mat` in the working directory. 

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
    * You should use `--selectedaks <n>` to choose the number of coefs (default: 5%)
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
    * **WARNING**: Make sure that the RPM level in `Mathcard.mo` is correct, and that time discretization in 
      `utils_openmodelica.runTestSimulation` matches the values used for the original dataset generation !!

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
		│   └── outputs/
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

1. Change paths in the script [Post-processing_Code/setupproj.m](Post-processing_Code/setupproj.m)
    * `output_path` must match the output path created by [test_driver.py](Deep_learning/test_driver.py)
    * `test_file_exact` must match any file in the path defined above
    
1. Run [Post-processing_Code/test_dnnmodelevaluation.m](Post-processing_Code/test_dnnmodelevaluation.m)
    * This creates `table.csv`, `Xexact.csv`, `Xpredicted.csv` files
    * `mkdir results` in the `output_folder` that contains all the data and move those files in the result folder
    
1. Run the notebook [notebooks/N7-Results-report.ipynb](notebooks/N7-Results-report.ipynb)
    * Setup the paths defined in the 2nd cell
    * This will create plots and format the result table
    * You can repeat all steps and set `SIMULATION_LVAD = False` in [Simulation_script/setup.py](Simulation_script/setup.py)
      to have a dataset without an assist device, the notebook will make plots to compare the two.

### Generate several datasets

The workflow described above works when we need to generate a single data set.
However, generating several datasets sequentially is pretty time consuming.
The solution is to run time-consuming steps in parallel.

For instance, step 3 is the most time-consuming step. Say that you want to generate 
datasets with different pump speed levels: 4000, 5000, 6000 RPM. You have to perform the
following:

a. Modify the pump speed in `Mathcard.mo` to 4000 RPM

b. Run steps 1-3

c. Repeat (a) and (b), but set RPM to 5000 and 6000

But one has to **be very careful once step 13 is reached**: the generation of the DNN test dataset
will be based on the current `Mathcard.mo` file, in which the RPM level has the value that was set
in (c) at last. 

To overcome this issue, at step 13, do not forget to re-modify `Mathcard.mo` with the RPM level that corresponds
to the current data set for which you are generating the DNN test dataset.

## Automatic pipeline

See documentation in `pipeline/script_pipeline.py`. 


