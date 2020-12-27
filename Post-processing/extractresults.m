function [variable, time] = extractresults(variableName,results)
% This function find a given parameter (VariableName) in the data file (results), 
% localises it, and returns its value (variable) in function of time (time).
%
% Input:
% VariableName: string. The name of the variable, e.g. LeftVentricle.PV
% results: Results from the simulation. This is a .mat file. 
% Description is in the README file.
%
% Output:
% variable: value of the variable
% time: time for each value.

%Define the matrix with the name, readable for human eye
results_name=results.name';
for i = 1:size(results_name,1)
    curnames = results_name(i,:);
    % search real caracters in current line
    findnonspaces = isstrprop(results_name(i,:),'alphanum') | ...
                    isstrprop(results_name(i,:),'punct');
    curnames = curnames(findnonspaces);
    if (strcmp(variableName,curnames))
        variableId = i;
    end
    
end

%Find informations in dataInfo. Takes only into account the first column.
variableInfo=results.dataInfo(:,variableId);
%Find if the data is in data set 1 or 2
variableDataSet=variableInfo(1);
%Find the position of the variable
variablePosition=abs(variableInfo(2));

%Time is the first row of data_2
variable=zeros(size(results.data_1,1),1);
%Extract the data in the right file (data_1 or data_2)
if variableDataSet==1
    variable=results.data_1(variablePosition,:);
    time=results.data_1(1,:);
    
elseif variableDataSet==2
    time=results.data_2(1,:);
    variable=results.data_2(variablePosition,:);
end
end