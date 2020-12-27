function [MeanArterialPressure] = MeanPressure(pressure)
% Compute the mean arterial pressure following.
% Mean arterial pressure (MAP) depends on systolic pressure (SAP) and 
% diastolic arterial pressure (DAP).
% MAP = DAP + 1/3*(SAP-DAP).

MeanArterialPressure=min(pressure)+1/3*(max(pressure)-min(pressure));
end