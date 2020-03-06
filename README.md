# VAD-0D-DNN
Internship with Jean Bonnemain & Simone Deparis

## Steps

1. Update the modelica code for Heart Mate III instead of Heart mate II. That is, fit polynomials for HQ curves of HMIII and modify the model:
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
