clear

% =========== PATH & FILE MANAGEMENT
setupproj

% Load the path where the data is. In this case, we want to load the results
% of the 0-D simulations.
FilePathOutputs = output_path; % from setupproj

filesExact = dir([FilePathOutputs,'*exact.mat']);
filesPredicted = dir([FilePathOutputs,'*predicted.mat']);

filesExact = sort_nat({filesExact.name});
filesPredicted = sort_nat({filesPredicted.name});

% For sensitivity analysis:
%filesExact = dir([FilePathOutputs,'*.mat']);
%filesPredicted = dir([FilePathOutputs,'*.mat']);


% =========== FIND SIMULATION TIME & TIME STEP
% Load first file to identify the time step and duration of the simulation
FirstFileName=test_file_exact; % from setupproj
res=load(FirstFileName);

TotalSimulationTime = (res.data_2(1,end));
tsub_min=ceil((1/2)*TotalSimulationTime);
tsub_max=TotalSimulationTime;


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
