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

MinExact=min(abs(Xexact));
MaxExact=max(abs(Xexact));
MinPredicted=min(abs(Xpredicted));
MaxPredicted=max(abs(Xpredicted));

MeanExact = mean(Xexact);
MeanPredicted = mean(Xpredicted);
SDexact=std(Xexact);
SDpredicted=std(Xpredicted);
meanErrorMatrix=mean(abs(Xexact-Xpredicted));
meanRelErrorMatrix=mean(abs(Xexact-Xpredicted)./(abs(Xexact)+1e-2));
SDErrorMatrix=std(abs(Xexact-Xpredicted));
[CImin,CImax]=confidenceinterval(abs(Xexact-Xpredicted),0.95);

%[HR,SAPM,SAPS,SAPD,PAPM,PAPS,PAPD,LVEF, LVEDV, LVESV, CI]

table = [
    MinExact;
    MaxExact;
    MeanExact;
    SDexact;
    MinPredicted;
    MaxPredicted;
    MeanPredicted;
    SDpredicted;
    meanErrorMatrix;
    meanRelErrorMatrix;
    SDErrorMatrix; %sigma
    CImin;
    CImax;
];
writematrix(table,'table.csv');
writematrix(Xexact,'Xexact.csv');
writematrix(Xpredicted,'Xpredicted.csv');

%x=(abs(Xexact(:,7)-Xpredicted(:,7)));
%pd=fitdist(x,'normal')
%ci=paramci(pd)