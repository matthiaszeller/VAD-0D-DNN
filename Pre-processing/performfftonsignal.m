function [y,aks,bks,T,t,signal] = performfftonsignal(signal,t,dt)

% for some reason, even when a constant timestep is specified in
% openmodelica, some timesteps are put sometimes in the middle. Here we get
% rid of them
[t,id] = getconstanttimesteps(t,dt);
signal = signal(id);

% we get rid of duplicate times if any
[t,it] = unique(t);
signal = signal(it);

% compute T, discard last values (for some reason this must be done..)
tmin = t(1);
tmax = t(end);
T = tmax-tmin;

t = t(1:end-1);
signal = signal(1:end-1);

% perform fft
[y,aks,bks] = myfft(signal);