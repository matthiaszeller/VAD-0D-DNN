

# Pre-processing

## About

Pre-process simulation output files by extracting systemic arterial pressure and pulmonary arterial pressure, converting 
those signals to their frequency-domain representation and storing `X.mat` and `Y.mat` files. 

This folder contains both MATLAB and Python code. 
Code was initially in MATLAB, but has been transcribed in Python because of 

1. Sampling issues. The Openmodelica output does not generate data points that are equally spaced 
in time. MATLAB code is known to work for T=20s, dt = 0.04s. Python code is known to work for 
   T=30s, dt=0.015s.
   
1. Automation of pre-processing procedure in the by the script `pipeline/script_pipeline.py`

## Getting started

Setup preprocessing by modifying the file `setup_preprocessing.py` (variables to extract, sampling, path of simulation 
files).