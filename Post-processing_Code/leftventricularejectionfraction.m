function [LVEF] = LeftVentricularEjectionFraction(LVVolume)
% Compute the left ventricular ejection fraction (LVEF).
% EDV: end diastolic ventricular volume
% ESV: end systolic ventricular volume
% LVEF = (EDV-ESV)/EDV

LVEF=((max(LVVolume)-min(LVVolume))/max(LVVolume))*100;
end