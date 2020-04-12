function [HR,period] = heartrate(data,time)
%Compute Heart rate and period

[pks,locs] = findpeaks(data);
period = max(diff(time(locs)));
HR=60*(1/period);
end
