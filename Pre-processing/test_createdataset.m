clear all
clc

% =========================== ABOUT
% Run this script once simulation data has been generated. 
% This generates the dataset for the deep learning network.
% Predictors: Fourier coefficients of the pulmonary artery pressure (PAP) 
% and systemic artery pressure (PAS). 
% Responses: cardiovascular parameters reflecting heart disfunction

% =========================== SETUP
% Setup the paths and time discretization
setup

% =========================== PROCESS PREDICTORS

% create input matrix
inputvariablesname = {'SystemicArteries.PC'; ...
                      'PulmonaryArteries.PC'};
nvariables = size(inputvariablesname,1);

files = dir([pathmats,'*.mat']);
nfiles = size(files,1);
% Important: reorder the output files in the order of integers
% E.g. {'file1', 'file13, 'file2'} becomes {'file1', 'file2', 'file13'}
% This is critical since the responses are ordered in the natural order
files = sort_nat({files.name});
count = 1;
for file = files
    disp(['Parsing input file ', num2str(count),'/',num2str(nfiles), ...
    ' ', file{1}]);
    results = load([pathmats,file{1}]);
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
    
    count = count + 1;
end

save('X.mat','X')

% =========================== PROCESS RESPONSES

disp('Extracting responses...')
% Since the predictors are sorted in the natural order of the samples,
% we can simply convert the parameters text files into a matrix file
Y = csvread(pathparams);
% Get rid of:
% - the first column 'n'
% - the first line, i.e. the header (column names)
Y = Y(2:end, 2:end);
save('Y.mat', 'Y')

disp('Done...')




