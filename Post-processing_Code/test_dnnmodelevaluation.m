clear

% =========== PATH & FILE MANAGEMENT
setupproj
% Load the path where the data is. In this case, we want to load the results
% of the 0-D simulations.
%FilePathOutputs='/Users/jean.bonnemain/Documents/Code/0d_model/Deep_learning/2020_01_13/outputs/';
FilePathOutputs = output_path; % from setupproj

% We also want to load the parameters.
FilePathParameters='/Users/jean.bonnemain/Documents/Code/0d_model/Deep_learning/2020_01_13/models/';

filesExact = dir([FilePathOutputs,'*exact.mat']);
filesPredicted = dir([FilePathOutputs,'*predicted.mat']);


% =========== FIND SIMULATION TIME & TIME STEP
% Load first file to identify the time step and duration of the simulation
FirstFileName=test_file_exact; % from setupproj
res=load(FirstFileName);

TotalSimulationTime = (res.data_2(1,end));
tsub_min=ceil((1/2)*TotalSimulationTime);
tsub_max=TotalSimulationTime;

% =========== DEFINE VARIABLES
% Names of the variables we want to analyze. Names come from the 0D model.
variablesname = {'SystemicArteries.PC','SAP'; ...
                      'SystemicArteries.PC','SAP'; ...
                      'PulmonaryArteries.PC', 'PAP'; ...
                      'LeftVentricle.V', 'LVV'; ...
                      'LeftVentricle.PV', 'LVP'; ...
                      'RightVentricle.PV', 'RVP'; ...
                      'RightVentricle.V', 'RVV'; ...
                      'LeftAtrium.PC', 'LAP'; ...
                      'AorticValve.Q', 'AoQ'; ...
                      };
                  

% Values to extract and analyze :
% LVEDV: Left ventricular end diastolic volume
% LVESV: Left ventricular end systolic volume
% LVEF: Left ventricular ejection fraction
% Heart rate
% SAPM, SAPS, SAPD
% PAPM, PAPS, PAPD
% Aortic flow, cardiac index

[Xexact]=createdataforanalysis(filesExact,FilePathOutputs,tsub_min, tsub_max);
[Xpredicted]=createdataforanalysis(filesPredicted,FilePathOutputs,tsub_min, tsub_max);

meanErrorMatrix=mean(abs(Xexact-Xpredicted));
meanRelErrorMatrix=mean(abs(Xexact-Xpredicted)./(abs(Xexact)+1e-2));
SDexact=std(Xexact);
SDpredicted=std(Xpredicted);
SDErrorMatrix=std(abs(Xexact-Xpredicted));
table = [SDErrorMatrix; meanErrorMatrix; mean(Xexact); meanRelErrorMatrix; SDexact];



% nfilesExact = size(filesExact,1);
% %nvariables = size(variablesname,1);
% nvariables=11;
% count = 1;
% Xexact = zeros(nfilesExact,nvariables);


% for file = filesExact'
%     results = load([FilePathOutputs,file.name]);
% %     for i=1:size(variablesname,1)
% %         variablename=variablesname{i,1};
% %         [var,t]=extractresults(variablename,results);
% %         [varsubrange,tsubrange]=timerange(var,t,tsub_min,tsub_max);
% %         data(i,:)=varsubrange;
% %     end
%     
%     
%     [SAP,t]=extractresults('SystemicArteries.PC',results);
%     [SAPsubrange,tsubrange]=timerange(SAP,t,tsub_min,tsub_max);
%     
%     SAPM = meanpressure(SAPsubrange);
%     SAPS=max(SAPsubrange);
%     SAPD=min(SAPsubrange);
%     
%     [HR,period] = heartrate(SAPsubrange,tsubrange);
%     
%     [PAP,t]=extractresults('PulmonaryArteries.PC',results);
%     [PAPsubrange,tsubrange]=timerange(PAP,t,tsub_min,tsub_max);
% 
%     PAPM = meanpressure(PAPsubrange);
%     PAPS=max(PAPsubrange);
%     PAPD=min(PAPsubrange);
%     
%     [LVV,t]=extractresults('LeftVentricle.V',results);
%     [LVVsubrange,tsubrange]=timerange(LVV,t,tsub_min,tsub_max);
%     
%     LVEF=leftventricularejectionfraction(LVVsubrange);
%     LVEDV=max(LVVsubrange);
%     LVESV=min(LVVsubrange);
%     
%     [AoQ,t]=extractresults('AorticValve.Q',results);
%     [AoQsubrange,tsubrange]=timerange(AoQ,t,tsub_min,tsub_max);
% 
% % Compute mean aortic flow (in L/min) and cardiac index (normalized for 
% %body surface, here 1.8m^2)
%     MeanAoQ = mean(AoQsubrange)*60/1000;
%     CI = MeanAoQ/1.8;
% 
% %     [LVP,t]=extractresults('LeftVentricle.PV',results);
% %     [LVPsubrange,tsubrange]=timerange(LVP,t,tsub_min,tsub_max);
% %     
% %     [RVV,t]=extractresults('RightVentricle.V',results);
% %     [RVVsubrange,tsubrange]=timerange(RVV,t,tsub_min,tsub_max);
% %     
% %     [RVP,t]=extractresults('RightVentricle.PV',results);
% %     [RVPsubrange,tsubrange]=timerange(RVP,t,tsub_min,tsub_max);
%     
%     Xexact(count,:)=[HR,SAPM,SAPS,SAPD,PAPM,PAPS,PAPD,LVEF, LVEDV, LVESV, CI];
%     count = count + 1;
% end

% nfilesPredicted = size(filesPredicted,1);
% Xpredicted = zeros(nfilesPredicted,nvariables);
% count = 1;
% 
% for file = filesPredicted'
%     results = load([FilePathOutputs,file.name]);
%     [SAP,t]=extractresults('SystemicArteries.PC',results);
%     [SAPsubrange,tsubrange]=timerange(SAP,t,tsub_min,tsub_max);
%     SAPM = meanpressure(SAPsubrange);
%     Xpredicted(count,1)=SAPM;
%     count = count + 1;
% end



