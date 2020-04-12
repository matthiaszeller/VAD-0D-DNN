% Script for the MATLAB project. Please read the file
% README.docx for more information.

% Load results
%FileName='Ursino1998Model_SHF_20s.mat';
%FileName='/Users/jean.bonnemain/Documents/Code/0d_model/Deep_learning/2020_01_13/outputs/Ursino1998Model_output_0exact.mat';

setupproj
Results=load(test_file_exact);

% Find the total simulation time in data_2
TotalSimulationTime = (Results.data_2(1,end));

% Define the interval of time we want to use to analyse data
% This because the first seconds of the simulation correspond to the 
% initiation of he simulation and do not have a physiological meaning.
Tsub_min=ceil((2/3)*TotalSimulationTime);
Tsub_max=TotalSimulationTime;

% Use the function "extractresults" to load a specific parameter in data.

% Systemic arterial pressure
[SAP,t]=extractresults('SystemicArteries.PC',Results);
[SAPsubrange,tsubrange]=timerange(SAP,t,Tsub_min,Tsub_max);

% Extract systemic systolic (SAP) and diastolic (DAP) values.
SAP=max(SAPsubrange);
DAP=min(SAPsubrange);

% Compute mean systemic arterial pressure.
MAP = meanpressure(SAPsubrange);

% Compute heart rate and period from the arterial curve
[HR,period] = heartrate(SAPsubrange,tsubrange);

% Pulmonary arterial pressure
[PAP,t]=extractresults('PulmonaryArteries.PC',Results);
[PAPsubrange,tsubrange]=timerange(PAP,t,Tsub_min,Tsub_max);

% Extract pulmonary systolic and diastolic values.
SPAP=max(PAPsubrange);
DPAP=min(PAPsubrange);

% Compute mean pulmonary arterial pressure.
MPAP=meanpressure(PAPsubrange);

% For left and right ventricle pressure and volume,
% we extract data for only one cardiac cycle,
% in order to plot the PV curves below.

% Left ventricular volume
[V_LV,t]=extractresults('LeftVentricle.V',Results);
[V_LVsubrange,tsubrangePVCurve]=timerange(V_LV,t,Tsub_max-period,Tsub_max);

% Extract left ventricular end diastolic volume and end systolic volume
EDV=max(V_LVsubrange);
ESV=min(V_LVsubrange);

% Compute left ventricular ejection fraction
LVEF=leftventricularejectionfraction(V_LVsubrange);

% Left ventricular pressure
[P_LV,t]=extractresults('LeftVentricle.PV',Results);
[P_LVsubrange,tsubrangePVCurve]=timerange(P_LV,t,Tsub_max-period,Tsub_max);

% Right ventricular pressure
[P_RV,t]=extractresults('RightVentricle.PV',Results);
[P_RVsubrange,tsubrangePVCurve]=timerange(P_RV,t,Tsub_max-period,Tsub_max);

% Right ventricular volume
[V_RV,t]=extractresults('RightVentricle.V',Results);
[V_RVsubrange,tsubrangePVCurve]=timerange(V_RV,t,Tsub_max-period,Tsub_max);

% Left atrial pressure
[L_Atria,t]=extractresults('LeftAtrium.PC',Results);
[L_Atriasubrange,tsubrange]=timerange(L_Atria,t,Tsub_min,Tsub_max);

% Compute Pulmonary Capillary Wedge Pressure
PCPW=mean(L_Atriasubrange);

% Flow through aortic valve
[AoQ,t]=extractresults('AorticValve.Q',Results);
[AoQsubrange,tsubrange]=timerange(AoQ,t,Tsub_min,Tsub_max);

% Compute mean aortic flow (in L/min) and cardiac index (normalized for 
%body surface, here 1.8m^2)
MeanAoQ = mean(AoQsubrange)*60/1000;
CI = MeanAoQ/1.8;

% Plot and save graphs for qualitative analysis.

% Plot Pressure - Volume Curve for left and right ventricles, and save
% it in a file (PVCurves.eps).
figure;
hold on;
plot(V_LVsubrange,P_LVsubrange, '-r');
plot(V_RVsubrange,P_RVsubrange, '-b');
xlabel('Volume [ml]','FontSize',16)
ylabel('Pressure [mmHg]','FontSize',16)
title('\it{Pressure - Volume Curves}','FontSize',22)
legend('Left Ventricle','Right Ventricle', 'Location','northwest');
set(gcf,'PaperPositionMode','auto')
print('PVCurves.eps','-depsc','-r0')
hold off;

% Plot systemic arterial pressure (only 2 seconds), and save it in a file
% called SystemicArterialPressure.eps.
figure;
hold on;
plot(tsubrange,SAPsubrange,'-r');
xlabel('Time [s]','FontSize',16)
ylabel('Pressure [mmHg]','FontSize',16)
title('\it{Systemic arterial pressure}','FontSize',22)
xlim([TotalSimulationTime-2,TotalSimulationTime]);
ylim([0,max(SAPsubrange)+10]);
set(gcf,'PaperPositionMode','auto')
print('SystemicArterialPressure.eps','-depsc','-r0')
hold off;

% Pulmonary arterial pressure (only 2 seconds), and save it in a file
% called PulmonaryArterialPressure2.eps.
figure;
hold on;
plot(tsubrange,PAPsubrange,'-r');
xlabel('Time [s]','FontSize',16)
ylabel('Pressure [mmHg]','FontSize',16)
title('\it{Pulmonary arterial pressure}','FontSize',22)
xlim([TotalSimulationTime-2,TotalSimulationTime]);
ylim([0,max(PAPsubrange)+10]);
set(gcf,'PaperPositionMode','auto')
print('PulmonaryArterialPressure2.eps','-depsc','-r0')
hold off;

% Here we compare the values of the experimental model (numerical) 
% with references obtained from Cox et al., Artif Organs 2009;33:593–603.
% And write all of them in a .csv file

% Reference values
EF_ref='15-33';
HR_ref='76-103';
CI_ref='1.9-2.4';
EDV_ref='205-522';
ESV_ref='140-249';
SAP_ref='107-115';
DAP_ref='68-76';
MAP_ref='78-95';
SPAP_ref='54-62';
DPAP_ref='28-29';
MPAP_ref='27-40';
PCWP_ref='17-29';

ReferenceValues={...
'15-33';
'76-103';
'1.9-2.4';
'205-522';
'140-249';
'107-115';
'68-76';
'78-95';
'54-62';
'28-29';
'27-40';
'17-29';};

% Text
Parameters={...
'LVEF';
'Heart rate';
'Cardiac index';
'End diastolic volume';
'End systolic volume';
'Systolic arterial pressure';
'Diastolic arterial pressure';
'Mean arterial pressure';
'Systolic pulmonary arterial pressure';
'Diastolic pulmonary arterial pressure';
'Mean pulmonary arterial pressure';
'Pulmonary capillary wedge pressure';};

% Experimental values
ExperimentalValues=[...
LVEF;
HR;
CI;
EDV;
ESV;
SAP;
DAP;
MAP;
SPAP;
DPAP;
MPAP;
PCPW;
];

% Units
Units={...
'[%]';
'[1/s]';
'[L/min]';
'[ml]';
'[ml]';
'[mmHg]';
'[mmHg]';
'[mmHg]';
'[mmHg]';
'[mmHg]';
'[mmHg]';
'[mmHg]';};

% Experimental values are "rounded" to have a precision of 2 digits. 
ExperimentalValues=round(ExperimentalValues*100)/100;

% Write data in a csv file. 
T=table(Parameters,ExperimentalValues,ReferenceValues,Units);
writetable(T,'ResultsSimulation.csv','WriteRowNames',true);