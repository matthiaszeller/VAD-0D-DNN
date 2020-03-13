function [coefs,aks,bks] = myfft(signal)

% we follow the algorithm presented in Scientific Computing with MATLAB and
% octave - Quarteroni, Saleri, Gervasio (4th edition). ref formula 3.26

if (size(signal,1) > 1)
    signal = signal';
end

n = length(signal) - 1;
h = 2 * pi / (n+1);

mu = mod(n,2);
M = (n-mu) / 2;

coefs = zeros(n+1,1);

indices = 1:n+1;

for ind = indices
    newValue = (exp(-1i * (indices - 1) * (ind - 1 - M) * h) * signal') / (n + 1);
    coefs(ind+mu) = newValue;
end

% fix first and last value for odd number of points
if (mu == 1)
    newValue = 1 / (2*(n+1)) * (-1).^(indices-1) * signal';
    coefs(1) = newValue;
    coefs(end) = newValue;
end

% % recover coefficients aks and bks for the even case
if (mu == 0)
    midindex = ceil((n+1)/2);
    aks = zeros(M+1,1);
    bks = zeros(M+1,1);
    aks(1) = coefs(midindex) * 2;
    bks(1) = 0;
    for ind = 1:M
        aks(ind+1) = coefs(midindex-ind) + coefs(midindex+ind);
        bks(ind+1) = 1i * (coefs(midindex+ind) - coefs(midindex-ind));
    end
else
    midindex = ceil((n+2)/2);
    aks = zeros(M+2,1);
    bks = zeros(M+2,1);
    aks(1) = coefs(midindex) * 2;
    bks(1) = 0;
    for ind = 1:M
        aks(ind+1) = coefs(midindex-ind) + coefs(midindex+ind);
        bks(ind+1) = 1i * (coefs(midindex+ind) - coefs(midindex-ind));
    end
    aks(end) = 2 * coefs(end);
    bks(end) = 0;
end

end

% Notes Jean about myfft.m
% 
% Ligne 21: ndices-1 : doit partir de 0 (j=0). Ind-1-M : -1-M permet de 
% partir de -M (car indices comments à 1). Signal’ : pas besoin de faire 
% boucle sur j, produit scalaire déjà effectué.
% Ligne 32 : equation 3.22.