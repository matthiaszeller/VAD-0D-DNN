clear all
close all
clc

% Load results
results=load('model_file_0.mat');

% The total time of the results is known and defined during the simulation.
% For this case, it is 300 seconds.
total_time=300;

% We want to show for a few heartbeats, i.e. few seconds.
tsub_min=5;
tsub_max=total_time;


% Extract data from the file, and then limit them to a given time 
% interval, e.g., Tsub_max-Tsub_min

% Left ventricular pressure
[P_LV,t]=extractresults('LeftVentricle.PV',results);
[P_LVsubrange,tsubrange]=timerange(P_LV,t,tsub_min,tsub_max);

% Left ventricular volume
[V_LV,t]=extractresults('LeftVentricle.V',results);
[V_LVsubrange,tsubrange]=timerange(V_LV,t,tsub_min,tsub_max);

% Right ventricular pressure
[P_RV,t]=extractresults('RightVentricle.PV',results);
[P_RVsubrange,tsubrange]=timerange(P_RV,t,tsub_min,tsub_max);

% Right ventricular volume
[V_RV,t]=extractresults('RightVentricle.V',results);
[V_RVsubrange,tsubrange]=timerange(V_RV,t,tsub_min,tsub_max);

% Flow through aortic valve
%[AoQ,t]=ExtractResults('AorticValve.Q',Results);
%[AoQsubrange,tsubrange]=TimeRange(AoQ,t,Tsub_min,Tsub_max);

% Systemic arterial pressure
[SAP,t]=extractresults('SystemicArteries.PC',results);
[SAPsubrange,tsubrange]=timerange(SAP,t,tsub_min,tsub_max);

% Pulmonary arterial pressure
[PAP,t]=extractresults('PulmonaryArteries.PC',results);
[PAPsubrange,tsubrange]=timerange(PAP,t,tsub_min,tsub_max);

% Plot Pressure - Volume Curve for left and right ventricles
figure;
hold on;
plot(V_LVsubrange,P_LVsubrange, '-r');
plot(V_RVsubrange,P_RVsubrange, '-b');
xlabel('Volume [ml]','FontSize',16)
ylabel('Pressure [mmHg]','FontSize',16)
title('\it{Pressure - Volume Curves}','FontSize',22)
legend('Left Heart','Right Heart', 'Location','northwest');
hold off;

% Compute left ventricular ejection fraction
LVEF=((max(V_LVsubrange)-min(V_LVsubrange))/max(V_LVsubrange))*100;

% Extract systemic systolic and diastolic values.
SAP=max(SAPsubrange);
DAP=min(SAPsubrange);

% Compute mean systemic arterial pressure.
Aortic_P_mean=integralmean(SAPsubrange,tsubrange);
MAP=min(SAPsubrange)+1/3*(max(SAPsubrange)-min(SAPsubrange));

% Extract pulmonary systolic and diastolic values.
SPAP=max(PAPsubrange);
DPAP=min(PAPsubrange);

% Compute mean pulmonary arterial pressure.
Pulmonary_P_mean=integralmean(PAPsubrange,tsubrange);
MPP=min(PAPsubrange)+1/3*(max(PAPsubrange)-min(PAPsubrange));

% Plot systemic arterial pressure
figure;
hold on;
plot(tsubrange,SAPsubrange,'-r');
xlabel('Time [s]','FontSize',16)
ylabel('Pressure [mmHg]','FontSize',16)
title('\it{Systemic arterial pressure}','FontSize',22)
saveas(gcf,'SAP.jpeg')
hold off;

% Plot pulmonary arterial pressure
figure;
hold on;
plot(tsubrange,PAPsubrange,'-r');
xlabel('Time [s]','FontSize',16)
ylabel('Pressure [mmHg]','FontSize',16)
title('\it{Pulmonary arterial pressure}','FontSize',22)
hold off;