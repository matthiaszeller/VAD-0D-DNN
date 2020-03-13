function [t,id] = getconstanttimesteps(t,dt)
% We use this function applied to a discretized time range t with constant
% dt, where are some intermediate timesteps. We discard these intermediate
% timesteps

tcopy = t - t(1);
id = find(abs(tcopy/dt-round(tcopy/dt))<1e-12);
t = t(id);