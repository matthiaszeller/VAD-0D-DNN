clear all
close all
clc

results=load('model_file_1.mat');
tsub_min = 10;
tsub_max = 20;

% take desired signal from results
[signal,t]=extractresults('SystemicArteries.PC',results);
[signal,t]=timerange(signal,t,tsub_min,tsub_max);

[y,aks,bks,T,t,signal] = performfftonsignal(signal,t,0.04);

% we reconstruct the original signal
reconstructedSignal = t * 0 + aks(1)/2;
for k = 2:length(aks)
    reconstructedSignal = reconstructedSignal + aks(k) * cos((k-1) * t * (2 * pi) / T) ...
                                              + bks(k) * sin((k-1) * t * (2 * pi) / T);
end

plot(t,signal,'.-r','Linewidth',2,'Markersize',15);
hold on
plot(t,reconstructedSignal,'--b','Linewidth',2)
plot(t,interpft(signal,length(signal)))
xlim([t(1) t(end)])
xlabel('t')
ylabel('SystemicArteries.PC')
set(gca,'fontsize', 15);
axis square
legend('signal','reconstructed','Location','southeast')
