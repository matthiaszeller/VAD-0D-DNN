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
1. Neural network training 
1. Neural network testing, statistics 

## Plan

1. ###  HQ curve

Update the modelica code for Heart Mate III since it is currently implemented for Heart Mate II.
	
* Fit quadratic polynomials to the HQ curves in the HMIII manual
* Duplicate the `VAD` model to `VAD2` and modify the expression of `Q` as a function of `dP`
* Duplicate the main model that instanciates the `VAD` model and instanciate `VAD2` instead
* Check the consistency of curves for medium heart failure (MHF) / severe heart failure (SHF) 
	and in context of baroreceptor regulation / no regulation. 
	Absence of regulation is simulated by setting the following parameters to zero: 

## Data

#### Modelica

* [Mathcard.mo](modelica/Mathcard.mo): modelica file containing the suitable libraries to model the cardiovascular system

#### Notebooks

1. [HQ curve HMIII](notebooks/N1-HQ-curve-HMIII.ipynb): fit quadratic polynomials of the pressure head-flow curves for different RPMs of the Heart Mate III
