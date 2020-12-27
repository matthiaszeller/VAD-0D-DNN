# VAD-0D-DNN

## About

This repository contains the code for the paper *Deep neural network to accurately predict left ventricular systolic 
function under mechanical assistance*, Bonnemain *et al.*

## Requirements

* Matlab R2020a
* [Openmodelica](https://openmodelica.org/) 1.16.0
* Python 3
    * Standard packages: matplotlib, pandas, numpy
    * [OMPython](https://github.com/OpenModelica/OMPython)
    * [Keras](https://keras.io/)


## Project Overview

The project is subdivided in the following sections:

* 0D Model Implementation in `modelica/original/Mathcard.mo`: 0D Model of the cardiovascular system implemented in the 
Modelica language.

* 0D Model Simulation in `Simulation_script`: Python code running the simulations with `OMPython` package
    * Generate a dataset by solving the model with different parameter values
    * Perform a sensitivity analysis: fix 3 heart-failure parameters and vary the 4th one
    * Case study: simulate a few scenarios by varying heart failure severity, pump speed

* Preprocessing in `Pre-processing`
    * Extract systemic arterial pressure and pulmonary arterial pressure from simulation data
    * Compute Fourier coefficients
    * Generate DNN input and output files

* DNN Training and Testing in `Deep_learning`: Python and (Keras Library)
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
    * `script_pipeline.py`: Automatic way to generate datasets with different pump configurations, merge them, train different DNNs on the 
    merged dataset
    *  `script_check_sim_data.py`: if an error occured while generating simulation data with `script_pipeline.py`, this 
    script helps you to identify the corrupted parts and prevents from recomputing everyting
    * `script_test_dnn.py` : once the best DNN is chosen, use this script to run simulations on test data with this 
    chosen DNN

## Reproductibility

The procedure described below was tested on Ubuntu 20.04. The settings as set in `script_pipeline.py` generated about 
2.5 TB of data. 

### Overview

1. Set the range of LVAD rotational speeds, the simulation settings (e.g. sampling), the number of Fourier 
coefficients to retain, the DNN architectures (and many other parameters) in `pipeline/script_pipeline.py`

1. Set sampling settings in `Pre-processing/setup_preprocessing.py`

1. Run the script `pipeline/script_pipeline.py`, which
    1. Generates simulation data for each specified pump configuration
    1. Merges data of all configurations in a single `.mat` file and creates training and test sets
    1. Trains the specified DNNs
    
    Note 1: this script takes several hours to run
    
    Note 2: the output folder will have the following contents:
    ```
    .
    ├── 4600_LVAD_AP
    │   ├── outputs
    │   │   ├── Ursino1998Model_VAD2_output_0.mat
    │   │   ├── Ursino1998Model_VAD2_output_1.mat
    │   │   ├── ...
    │   ├── parameters.txt
    │   ├── X.mat
    │   ├── Y.mat
    │   └── ...
    ├── 4700_LVAD_AP
    │   ├── ...
    ├── ...
    ├── dataset.mat
    ├── dnns
    │   ├── dnn_10_layers_128_neurons
    │   │   ├── DNN_Performance.eps
    │   │   ├── history.bin
    │   │   ├── losses.eps
    │   │   ├── model.h5
    │   │   ├── Ytestpred.txt
    │   │   └── Ytest.txt
    │   ├── dnn_10_layers_16_neurons
    │   │   ├── ...
    │   ├── ...
    │   └── normdata.mat
    ├── logs.log
    └── mainlog.log
    ```
   
1. Once you identified a DNN architecture that performs well, 
run `pipeline/script_test_dnn.py` to generate simulations on test data

1. Compute hemodynamic quantities on test data:
    1. Setup paths in `Post-processing_Code/setupproj.m` to target the `outputs/` folder of the selected DNN architecture
    1. Run `test_dnnmodelevaluation.m`
