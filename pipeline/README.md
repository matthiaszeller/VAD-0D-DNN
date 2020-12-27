
# Pipeline

## About

Wraps simulation data generation, preprocessing, and DNN training in the script `script_pipeline.py`. 
All settings are defined in `script_pipeline.py` and `Pre-processing/setup_preprocessing.py`. 
Once you identified a DNN architecture that performs well, run `script_test_dnn.py` to generate simulation based on 
0D-model parameters predicted by the DNN, and use post-processing code (in `Post-processing`) to compute hemodynamic 
quantities and their errors on the test set. 