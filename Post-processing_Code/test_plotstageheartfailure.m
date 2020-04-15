clear
close all

setupproj

FilePathOutputs=output_path;
files = dir([FilePathOutputs,'*.mat']);
numfiles=length(files)
% Allocate an array to contain all the data from output files
mydata=cell(numfiles)

for i=1:numfiles
    mydata{i}=load([FilePathOutputs,files(i).name]);
end

% Find the total simulation time in data_2
TotalSimulationTime = (mydata{1}.data_2(1,end));

% Define the interval of time we want to use to analyse data
% This because the first seconds of the simulation correspond to the 
% initiation of he simulation and do not have a physiological meaning.
Tsub_min=ceil((2/3)*TotalSimulationTime);
Tsub_max=TotalSimulationTime;

% Use the function "extractresults" to load a specific parameter in data.

for i=1:numfiles
    disp(i)
    Results=mydata{i};
    % Systemic arterial pressure
    [data(i).SAP,t]=extractresults('SystemicArteries.PC',Results);
    [data(i).SAPsubrange,tsubrange]=timerange(data(i).SAP,t,Tsub_min,Tsub_max);

    % Extract systemic systolic (SAP) and diastolic (DAP) values.
    data(i).SAP=max(data(i).SAPsubrange);
    data(i).DAP=min(data(i).SAPsubrange);

    % Compute mean systemic arterial pressure.
    data(i).MAP = meanpressure(data(i).SAPsubrange);

    % Compute heart rate and period from the systemic arterial curve
    [data(i).HR,data(i).period] = heartrate(data(i).SAPsubrange,tsubrange);

    % Pulmonary arterial pressure
    [data(i).PAP,t]=extractresults('PulmonaryArteries.PC',Results);
    [data(i).PAPsubrange,tsubrange]=timerange(data(i).PAP,t,Tsub_min,Tsub_max);

    % Extract pulmonary systolic and diastolic values.
    data(i).SPAP=max(data(i).PAPsubrange);
    data(i).DPAP=min(data(i).PAPsubrange);

    % Compute mean pulmonary arterial pressure.
    data(i).MPAP=meanpressure(data(i).PAPsubrange);

    % For left and right ventricle pressure and volume,
    % we extract data for only one cardiac cycle,
    % in order to plot the PV curves below.

    % Left ventricular volume
    [data(i).V_LV,t]=extractresults('LeftVentricle.V',Results);
    [data(i).V_LVsubrange,tsubrangePVCurve]=timerange(data(i).V_LV,t,Tsub_max-data(i).period,Tsub_max);

    % Extract left ventricular end diastolic volume and end systolic volume
    data(i).EDV=max(data(i).V_LVsubrange);
    data(i).ESV=min(data(i).V_LVsubrange);

    % Compute left ventricular ejection fraction
    data(i).LVEF=leftventricularejectionfraction(data(i).V_LVsubrange);

    % Left ventricular pressure
    [data(i).P_LV,t]=extractresults('LeftVentricle.PV',Results);
    [data(i).P_LVsubrange,tsubrangePVCurve]=timerange(data(i).P_LV,t,Tsub_max-data(i).period,Tsub_max);

    % Right ventricular pressure
    [data(i).P_RV,t]=extractresults('RightVentricle.PV',Results);
    [data(i).P_RVsubrange,tsubrangePVCurve]=timerange(data(i).P_RV,t,Tsub_max-data(i).period,Tsub_max);

    % Right ventricular volume
    [data(i).V_RV,t]=extractresults('RightVentricle.V',Results);
    [data(i).V_RVsubrange,tsubrangePVCurve]=timerange(data(i).V_RV,t,Tsub_max-data(i).period,Tsub_max);

    % Left atrial pressure
    [data(i).L_Atria,t]=extractresults('LeftAtrium.PC',Results);
    [data(i).L_Atriasubrange,tsubrange]=timerange(data(i).L_Atria,t,Tsub_min,Tsub_max);

    % Compute Pulmonary Capillary Wedge Pressure
    data(i).PCPW=mean(data(i).L_Atriasubrange);

    % Flow through aortic valve
    [data(i).AoQ,t]=extractresults('AorticValve.Q',Results);
    [data(i).AoQsubrange,tsubrange]=timerange(data(i).AoQ,t,Tsub_min,Tsub_max);

    % Compute mean aortic flow (in L/min) and cardiac index (normalized for 
    %body surface, here 1.8m^2)
    data(i).MeanAoQ = mean(data(i).AoQsubrange)*60/1000;
    data(i).CI_mean = data(i).MeanAoQ/1.8;
    data(i).MeanAoQ = meanintegral(data(i).AoQsubrange,tsubrange)*60/1000;
    data(i).CI_fct  = data(i).MeanAoQ/1.8;
end

% % Plot and save graphs for qualitative analysis.
% 
% % Plot Pressure - Volume Curve for left and right ventricles, and save
% % it in a file (PVCurves.eps).
% figure;
% hold on;
% plot(V_LVsubrange,P_LVsubrange, '-r');
% plot(V_RVsubrange,P_RVsubrange, '-b');
% xlabel('Volume [ml]','FontSize',16)
% ylabel('Pressure [mmHg]','FontSize',16)
% title('\it{Pressure - Volume Curves}','FontSize',22)
% legend('Left Ventricle','Right Ventricle', 'Location','northwest');
% set(gcf,'PaperPositionMode','auto')
% print('PVCurves.eps','-depsc','-r0')
% hold off;
% 
% % 
figure;
hold on;
plot(data(2).V_LVsubrange,data(2).P_LVsubrange, '-b', 'LineWidth',2);
plot(data(1).V_LVsubrange,data(1).P_LVsubrange, '-g','LineWidth',2);
plot(data(3).V_LVsubrange,data(3).P_LVsubrange, '-r','LineWidth',2);
legend('Normal Heart','Moderate Heart Failure','Severe Heart Failure');
xlabel('Volume [ml]','FontSize',16)
ylabel('Pressure [mmHg]','FontSize',16)
title('\it{Pressure - Volume Curves}','FontSize',22)
print('PVCurves.eps','-depsc','-r0')
hold off;
% 

figure;
hold on;
plot(tsubrange,data(3).SAPsubrange,'-k', 'LineWidth',3);
%plot(tsubrangepredicted,SAPPredictedsubrange,'--b');
xlabel('Time [s]','FontSize',16)
ylabel('Pressure [mmHg]','FontSize',16)
title('\it{Systemic arterial pressure}','FontSize',22)
xlim([TotalSimulationTime-2,TotalSimulationTime]);
ylim([0,max(data(3).SAPsubrange)+10]);
hold off;

figure;
hold on;
plot(data(2).V_LVsubrange,data(2).P_LVsubrange, '-m', 'LineWidth',2);
xlabel('Volume [ml]','FontSize',16)
ylabel('Pressure [mmHg]','FontSize',16)
title('\it{Pressure - Volume Curve}','FontSize',22)
xlim([0,160]);
print('PVCurves_global_concept.eps','-depsc','-r0')
hold off;