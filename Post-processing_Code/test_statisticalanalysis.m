clear all
clc

setupproj
pathmats = output_path;
%pathmats = '/Users/jean.bonnemain/Documents/Results/0d_model/2019_12_05/outputs/';

% create input matrix with predicted data
inputvariablesname = {'SystemicArteries.PC'; ...
                      'PulmonaryArteries.PC'};
nvariables = size(inputvariablesname,1);
% specifics of the time discretization
tsub_min = 15;
tsub_max = 20;
dt = 0.04;

% create input matrix with exact data
files = dir([pathmats,'*exact.mat']);
nfiles = size(files,1);
count = 1;
for file = files'
    results = load([pathmats,file.name]); 
    for i = 1:nvariables
        [signal,t]=extractresults(inputvariablesname{i},results);
        [signal,t]=timerange(signal,t,tsub_min,tsub_max);
        
        [t,id] = getconstanttimesteps(t,dt);
        signal = signal(id);
        
        % we get rid of duplicate times if any
        [t,it] = unique(t);
        signal = signal(it);
        if (count == 1 && i == 1)
            exactsignals = zeros(nfiles,nvariables,size(signal,2));
            exacttimesteps = zeros(nfiles,nvariables,size(signal,2));
        end
        exactsignals(count,i,:) = signal;
        exacttimesteps(count,i,:) = t;
        [y,aks,bks,T,t,signal] = performfftonsignal(signal,t,dt);
        % some coefficients in bks are always zero            
        bks = bks(abs(bks) > 1e-15);
        if (count == 1 && i == 1)
            n = length(aks) + length(bks);
            Xexact = zeros(nfiles,nvariables,n);
        end
        Xexact(count,i,:) = [aks;bks];
    end
    
    disp(['Parsing input file ', num2str(count),'/',num2str(nfiles)]);
    count = count + 1;
end

save('Xexact.mat','Xexact')

files = dir([pathmats,'*predicted.mat']);
nfiles = size(files,1);
count = 1;
for file = files'
    results = load([pathmats,file.name]); 
    for i = 1:nvariables
        [signal,t]=extractresults(inputvariablesname{i},results);
        [signal,t]=timerange(signal,t,tsub_min,tsub_max);
        
        [t,id] = getconstanttimesteps(t,dt);
        signal = signal(id);
        
        % we get rid of duplicate times if any
        [t,it] = unique(t);
        signal = signal(it);
        
        [c,lags] = xcorr(signal, squeeze(exactsignals(count,i,:)));
        
        [y,aks,bks,T,t,signal] = performfftonsignal(signal,t,dt);
        % some coefficients in bks are always zero            
        bks = bks(abs(bks) > 1e-15);
        if (count == 1 && i == 1)
            n = length(aks) + length(bks);
            Xpredicted = zeros(nfiles,nvariables,n);
        end
        Xpredicted(count,i,:) = [aks;bks];
    end
    
    disp(['Parsing input file ', num2str(count),'/',num2str(nfiles)]);
    count = count + 1;
end

save('Xpredicted.mat','Xpredicted')


Xdiff=Xpredicted-Xexact;
Xdiff(:,:,1)=0.5*Xdiff(:,:,1); % First coefficient is alpha/2
Xpredicted(:,:,1)=0.5*Xpredicted(:,:,1); % First coefficient is alpha/2
Errors=zeros(size(Xdiff,1),size(Xdiff,2));
Norms=zeros(size(Xdiff,1),size(Xdiff,2));
% Compute the norm for the differences between predicted and exact
% coefficients.
% Compute de L2 norm on the Fourier coefficients.
for i=1:size(Xdiff,1)
    for j=1:size(Xdiff,2) %for each parameter (inputvariablesname)
        Errors(i,j)=norm(squeeze(Xdiff(i,j,:)));% Norm for all coefficients of 1 simulation
        Norms(i,j)=norm(squeeze(Xpredicted(i,j,:))); % Norm for a sample
    end
end

Errors = Errors * sqrt(T); %Multiply by sqrt(T) (sqrt(N)*int from 0 to N (f(x)))
Norms = Norms * sqrt(T);
ErrorsTot = Errors ./ Norms; % Relative error
ErrorNormalizedParam1=mean(ErrorsTot(:,1));
ErrorNormalizedParam2=mean(ErrorsTot(:,2));


ArteryTest=load('/Users/jean.bonnemain/Documents/Results/0d_model/2019_12_05/outputs/Ursino1998Model_output_5exact.mat');
[SAPExact,t]=extractresults('SystemicArteries.PC',ArteryTest);
[SAPExactsubrange,tsubrangexact]=timerange(SAPExact,t,tsub_min,tsub_max);

figure;
hold on;
plot(tsubrangexact,SAPExactsubrange,'-r');
%plot(tsubrangepredicted,SAPPredictedsubrange,'--b');
xlabel('Time [s]','FontSize',16)
ylabel('Pressure [mmHg]','FontSize',16)
title('\it{Systemic arterial pressure}','FontSize',22)
%xlim([TotalSimulationTime-2,TotalSimulationTime]);
ylim([0,max(SAPExactsubrange)+10]);
hold off;

T= tsub_max-tsub_min;
Integration=sqrt(trapz(tsubrangexact,SAPExactsubrange.^2));
%Compare L2 norm on Fourier coef and using the trapz on the native
%function. Fourier coef are better.
% Pay attention to phase! Signals are almost the same but shifted. Use
% xcorr: doesn't work.
% Try cshift, calcule avec mesure des peaks à tsubmin

