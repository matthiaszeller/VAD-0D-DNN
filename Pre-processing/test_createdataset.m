clear all
clc

%Modify both paths for X and Y!

pathmats = '/Volumes/BONNEMAIN/2019_12_24/outputs/';

% create input matrix
inputvariablesname = {'SystemicArteries.PC'; ...
                      'PulmonaryArteries.PC'};
nvariables = size(inputvariablesname,1);
% specifics of the time discretization
tsub_min = 10;
tsub_max = 20;
dt = 0.04;

files = dir([pathmats,'*.mat']);
nfiles = size(files,1);
count = 1;
for file = files'
    results = load([pathmats,file.name]); 
    for i = 1:nvariables
        [signal,t]=extractresults(inputvariablesname{i},results);
        [signal,t]=timerange(signal,t,tsub_min,tsub_max);
        [y,aks,bks,T,t,signal] = performfftonsignal(signal,t,dt);
        % some coefficients in bks are always zero            
        bks = bks(abs(bks) > 1e-15);
        if (count == 1 && i == 1)
            n = length(aks) + length(bks);
            X = zeros(nfiles,nvariables,n);
        end
        X(count,i,:) = [aks;bks];
    end
    
    disp(['Parsing input file ', num2str(count),'/',num2str(nfiles)]);
    count = count + 1;
end

save('X.mat','X')

% create output matrix
pathmos = '/Volumes/BONNEMAIN/2019_12_24/models/';
model = 'ModelParametersNH';

listparameters = {'Param_LeftVentricle_Emax0'; ...
                  'Param_LeftVentricle_EmaxRef0'; ...
                  'Param_LeftVentricle_AGain_Emax'; ...
                  'Param_LeftVentricle_kE'};

files = dir([pathmos,'*.mo']);
nfiles = size(files,1);
Y = [];
count = 1;
for file = files'
    params = getparametersfrommodel([pathmos,file.name],listparameters,model);
    if (count == 1)
        n = length(params);
        Y = zeros(nfiles,n);
    end
    
    Y(count,:) = params;
  
    disp(['Parsing output file ', num2str(count),'/',num2str(nfiles)]);
    count = count + 1;
end

save('Y.mat','Y')