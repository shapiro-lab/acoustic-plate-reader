function [ guidelines ] = my_qc(EC50, Hills_Slope)
%Sigmoidal Plot QC
% From "Guidelines for accurate EC50/IC50 estimation", Sebaugh, JL; Pharmaceutical Statistics 2011
    % A guideline of "at least two concentrations are required beyond the lower
    % and upper bend points in order to accurately estimate the relative EC/IC50."

    % The beginning and end of the linear portion of the beginning (upper) and
    % end (lower) bend points of the linear portion of the curve are given by
    % upper  = EC50 * k^(1/Hills_Slope)
    % lower = EC50 * (1/k)^(1/Hills_Slope)
    % k = 4.6805

%%
% let's say 1 <---> 1 uM so we're ranging from 1 nM to 1 M, for example
% we could change X or even pass X into the function as needed
X = log10([0.001, 0.01, 0.1, 1, 10, 100, 1000, 10000, 100000, 1000000]);  

k = 4.6805;

lower_bend = EC50 * ((1/k)^(1/Hills_Slope));
upper_bend = EC50 * (k^(1/Hills_Slope));

lower_bend_limit = log10(lower_bend);
lower_bend_count = nnz(X < lower_bend_limit);

upper_bend_limit = log10(upper_bend);
upper_bend_count = nnz(X > upper_bend_limit);

guidelines = table(lower_bend, upper_bend, (lower_bend_count > 1.9)  &  (upper_bend_count > 1.9));
guidelines.Properties.VariableNames = {'lower_bend', 'upper_bend', 'Reliable'};
