function [X] = createdataforanalysis(files,path,tsub_min, tsub_max)
%CREATEDATAFORANALYSIS Load data from the 0D model to perform statistical
%analysis. Returns a vector where each row is a simulation, and columns are
%as follows: HR,SAPM,SAPS,SAPD,PAPM,PAPS,PAPD,LVEF, LVEDV, LVESV, CI.
%HR: Heart rate
%SAPM, S, D: Systemic Arterial Pressure Mean/Systoloc/Diastolic
%PAPM, S, D: Pulmonary Arterial Pressure Mean/Systoloc/Diastolic
%LVEF: Left Ventricular Ejection Fraction
%LVEDV: Left Ventricular End Diastolic Volume
%LVESV: Left Ventricular End Systolic Volume
%CI: Cardiac Index

%   Detailed explanation goes here
count=1;

for file = files'
    results = load([path,file.name]);
    
    [SAP,t]=extractresults('SystemicArteries.PC',results);
    [SAPsubrange,tsubrange]=timerange(SAP,t,tsub_min,tsub_max);
    
    SAPM = meanpressure(SAPsubrange);
    SAPS=max(SAPsubrange);
    SAPD=min(SAPsubrange);
    
    [HR,period] = heartrate(SAPsubrange,tsubrange);
    
    [PAP,t]=extractresults('PulmonaryArteries.PC',results);
    [PAPsubrange,tsubrange]=timerange(PAP,t,tsub_min,tsub_max);

    PAPM = meanpressure(PAPsubrange);
    PAPS=max(PAPsubrange);
    PAPD=min(PAPsubrange);
    
    [LVV,t]=extractresults('LeftVentricle.V',results);
    [LVVsubrange,tsubrange]=timerange(LVV,t,tsub_min,tsub_max);
    
    LVEF=leftventricularejectionfraction(LVVsubrange);
    LVEDV=max(LVVsubrange);
    LVESV=min(LVVsubrange);
    
    [AoQ,t]=extractresults('AorticValve.Q',results);
    [AoQsubrange,tsubrange]=timerange(AoQ,t,tsub_min,tsub_max);

    % Compute mean aortic flow (in L/min) and cardiac index (normalized for 
    %body surface, here 1.8m^2)
    MeanAoQ = meanintegral(AoQsubrange,tsubrange)*60/1000;
    CI = MeanAoQ/1.8;

    % Left atrial pressure
    [L_Atria,t]=extractresults('LeftAtrium.PC',results);
    [L_Atriasubrange,tsubrange]=timerange(L_Atria,t,tsub_min,tsub_max);

    % Compute Pulmonary Capillary Wedge Pressure
    PCPW=mean(L_Atriasubrange);
    
%     [LVP,t]=extractresults('LeftVentricle.PV',results);
%     [LVPsubrange,tsubrange]=timerange(LVP,t,tsub_min,tsub_max);
%     
%     [RVV,t]=extractresults('RightVentricle.V',results);
%     [RVVsubrange,tsubrange]=timerange(RVV,t,tsub_min,tsub_max);
%     
%     [RVP,t]=extractresults('RightVentricle.PV',results);
%     [RVPsubrange,tsubrange]=timerange(RVP,t,tsub_min,tsub_max);
    
    X(count,:)=[HR,SAPM,SAPS,SAPD,PAPM,PAPS,PAPD,LVEF, LVEDV, LVESV, CI, PCPW];
    count = count + 1;
end
end