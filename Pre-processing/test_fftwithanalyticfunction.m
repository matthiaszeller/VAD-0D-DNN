clear all
close all
clc

% define function to use as test
n = 26;
T = 2.4 * pi;
t = T / (n+1) * (0:n);
signal = t.*(t-2*pi).*exp(-t);

% perform fft
[y,aks,bks] = myfft(signal);

% we get the same result using this combination of matlab native functions
% y = conj(fftshift(fft(signal)))'/(n+1);

% reconstruct signal by following formula 3.24 in Scientific Computing 
% with MATLAB and octave - Quarteroni, Saleri, Gervasio (4th edition)
mu = mod(n,2);
M = floor((n-mu)/2);
reconstructedSignal = t * 0;
for k = -M:1:M
    reconstructedSignal = reconstructedSignal + y(k + (M+mu) + 1) * ...
                          exp(1i * k * t * (2 * pi) / T);
end
reconstructedSignal = reconstructedSignal + 2 * mu * y(end) * ...
                      cos((M+1) * t * (2 * pi) / T);

subplot(1,2,1)
% plot original, reconstructed and reconstructed signal 
% using matlab's interpft function 
plot(t,signal,'.-r','Linewidth',2,'Markersize',15);
hold on
plot(t,reconstructedSignal,'--b','Linewidth',2)
plot(t,interpft(signal,n+1),'--g','Linewidth',2);
xlim([0 t(end)])
xlabel('t')
ylabel('f(t)')
set(gca,'fontsize', 15);
legend('signal','reconstructed','reconstructed by matlab','Location','southeast')
axis square

% here we reconstruct using the coefficients aks and bks
reconstructedSignal = t * 0 + aks(1)/2;
for k = 2:length(aks)
    reconstructedSignal = reconstructedSignal + aks(k) * cos((k-1) * t * (2 * pi) / T) ...
                                              + bks(k) * sin((k-1) * t * (2 * pi) / T);
end

subplot(1,2,2)
% plot original, reconstructed and reconstructed signal 
% using matlab's interpft function 
plot(t,signal,'.-r','Linewidth',2,'Markersize',15);
hold on
plot(t,reconstructedSignal,'--b','Linewidth',2)
plot(t,interpft(signal,n+1),'--g','Linewidth',2);
xlim([0 t(end)])
xlabel('t')
ylabel('f(t)')
set(gca,'fontsize', 15);
legend('signal','reconstructed','reconstructed by matlab','Location','southeast')
axis square