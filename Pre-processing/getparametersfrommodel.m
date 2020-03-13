function [params] = getparametersfrommodel(filename,listparameters,model)

fid = fopen(filename,'r');
output = fgets(fid);
nparameters = size(listparameters,1);
searchforparams = -1;
params = zeros(nparameters,1);
while output ~= -1
    if contains(output,['model ',model])
        searchforparams = 0;
    end
    if (searchforparams >= 0 && searchforparams < nparameters)
        for ind = 1:nparameters
            if contains(output,listparameters{ind})
                lastequal = find(output == '=', 1, 'last');
                output = output(lastequal+2:end);
                lastesemicolon = find(output == ';', 1, 'last');
                output = output(1:lastesemicolon-1);
                searchforparams = searchforparams + 1;
                params(ind) = str2double(output);
            end
        end
    end
    if (searchforparams == nparameters)
        break;
    end
    output = fgets(fid);
end

fclose(fid);