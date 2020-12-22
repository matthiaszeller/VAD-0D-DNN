
# Deep learning

## About

Train deep neural networks to predict (left) heart-failure parameters from the truncated frequency-domain representation 
of systemic arterial pressure and pulmonary arterial pressure under left ventricular assist device (Heart Mate III) with 
Artificial Pulse. 


## DNN architecture selection

First, use the `script_stats_dnn.py` file to train several DNNs with different achitectures (number of layers, number of 
neurons per layer) as well as different input sizes (number of Fourier coefficients to retain). The resulting DNNs can 
be compared through the validation mean absolute error (to identify the best-performing DNN) and through the learning 
curves (to identify overfitting).

Use the notebook `N16-Analyze-DNN-stats.ipynb` to analyze the results. 

Note: it might be preferred to first select a subset of input sizes with the notebook 
`N15-Fourier-coefs-selection.ipynb`.

