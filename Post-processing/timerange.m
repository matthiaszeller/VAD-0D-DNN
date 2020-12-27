function [variableRange, time] = timerange(variable,t,t_min,t_max)
% This function extract the variable (Variable), for a given interval of
% time (t_min, t_max). It returns the value of the variable for this
% interval of time.
indices = find(t >= t_min & t <= t_max);
variableRange=variable(indices);
time=t(indices);
end