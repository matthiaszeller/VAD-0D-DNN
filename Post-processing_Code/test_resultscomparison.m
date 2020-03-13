clear all
close all
clc

FilePath='/Users/jean.bonnemain/Documents/Results/0d_model/2019_12_05/outputs/';
%FilePath='/Users/jean.bonnemain/Documents/Code/0d_model/Simulation_script/2019_12_05/outputs/';

resultsExact=load([FilePath,'Ursino1998Model_output_250exact.mat']);
resultsPredicted=load([FilePath,'Ursino1998Model_output_250predicted.mat']);

TotalSimulationTime = (resultsExact.data_2(1,end));
%tsub_min=ceil((2/3)*TotalSimulationTime);
tsub_min=0;
tsub_max=TotalSimulationTime;

[SAPExact,t]=extractresults('SystemicArteries.PC',resultsExact);
[SAPExactsubrange,tsubrangexact]=timerange(SAPExact,t,tsub_min,tsub_max);

[SAPPredicted,t]=extractresults('SystemicArteries.PC',resultsPredicted);
[SAPPredictedsubrange,tsubrangepredicted]=timerange(SAPPredicted,t,tsub_min,tsub_max);

[HRExact,periodExact] = HeartRate(SAPExactsubrange,tsubrangexact);
[HRPredicted,periodPredicted] = HeartRate(SAPPredictedsubrange,tsubrangepredicted);

% Left ventricular volume
[V_LVExact,t]=extractresults('LeftVentricle.V',resultsExact);
[V_LVExactsubrange,tsubrangePVCurve]=TimeRange(V_LVExact,t,tsub_max-periodExact,tsub_max);

% Left ventricular volume
[V_LVPredicted,t]=extractresults('LeftVentricle.V',resultsPredicted);
[V_LVPredictedsubrange,tsubrangePVCurve]=TimeRange(V_LVPredicted,t,tsub_max-periodPredicted,tsub_max);

% Left ventricular volume
[P_LVExact,t]=extractresults('LeftVentricle.PV',resultsExact);
[P_LVExactsubrange,tsubrangePVCurve]=TimeRange(P_LVExact,t,tsub_max-periodExact,tsub_max);

% Left ventricular volume
[P_LVPredicted,t]=extractresults('LeftVentricle.PV',resultsPredicted);
[P_LVPredictedsubrange,tsubrangePVCurve]=TimeRange(P_LVPredicted,t,tsub_max-periodPredicted,tsub_max);


% Extract end diastolic volume and end systolic volume
%EDV=max(V_LVsubrange);
%ESV=min(V_LVsubrange);

% Compute left ventricular ejection fraction
%LVEF=LeftVentricularEjectionFraction(V_LVsubrange);

% Flow through aortic valve
[AoQExact,t]=ExtractResults('AorticValve.Q',resultsExact);
[AoQExactsubrange,tsubrangexact]=TimeRange(AoQExact,t,tsub_min,tsub_max);

[AoQPredicted,t]=ExtractResults('AorticValve.Q',resultsPredicted);
[AoQPredictedsubrange,tsubrangepredicted]=TimeRange(AoQPredicted,t,tsub_min,tsub_max);


[SAPExact,t]=extractresults('SystemicArteries.PC',resultsExact);
[SAPExactsubrange,tsubrangexact]=timerange(SAPExact,t,tsub_min,tsub_max);

[SAPPredicted,t]=extractresults('SystemicArteries.PC',resultsPredicted);
[SAPPredictedsubrange,tsubrangepredicted]=timerange(SAPPredicted,t,tsub_min,tsub_max);


figure;
hold on;
plot(tsubrangexact,SAPExactsubrange,'-r');
plot(tsubrangepredicted,SAPPredictedsubrange,'--b');
xlabel('Time [s]','FontSize',16)
ylabel('Pressure [mmHg]','FontSize',16)
title('\it{Systemic arterial pressure}','FontSize',22)
xlim([TotalSimulationTime-2,TotalSimulationTime]);
ylim([0,max(SAPExactsubrange)+10]);
hold off;

figure;
hold on;
plot(tsubrangexact,AoQExactsubrange,'-r');
plot(tsubrangepredicted,AoQPredictedsubrange,'--b');
xlabel('Time [s]','FontSize',16)
ylabel('Flow [L/min]','FontSize',16)
title('\it{Aortic flow}','FontSize',22)
xlim([0,TotalSimulationTime]);
%ylim([0,max(AoQExactsubrange)]);
hold off;


% Plot Pressure - Volume Curve for left and right ventricles
figure;
hold on;
plot(V_LVExactsubrange,P_LVExactsubrange, '-r');
plot(V_LVPredictedsubrange,P_LVPredictedsubrange, '--b');
xlabel('Volume [ml]','FontSize',16)
ylabel('Pressure [mmHg]','FontSize',16)
title('\it{Pressure - Volume Curves}','FontSize',22)
%legend('Left Ventricle','Right Ventricle', 'Location','northwest');
hold off;