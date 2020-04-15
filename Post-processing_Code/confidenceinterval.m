function [CImin,CImax] = confidenceinterval(X,perc)
%CONFIDENCEINTERVAL Compute confidence interval. Z values: 80% - 1.28, 90%
%- 1.645, 95% - 1.96, 98% - 2.33, 99% - 2.58

if perc==0.8
    z=1.28;
elseif perc==0.9
        z=1.645;
elseif perc==0.95
            z=1.96;
elseif perc==0.98
                z=2.33;
else perc==0.99
                    z=2.58;
end

    CImin=mean(abs(X))-z*(std(X)/sqrt(size(X,1)))
    CImax=mean(abs(X))+z*(std(X)/sqrt(size(X,1)))
end

