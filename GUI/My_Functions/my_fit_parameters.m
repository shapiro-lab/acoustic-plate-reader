function [ fit_parameters_row ] = my_fit_parameters( fitresult, gof )
%UNTITLED12 Summary of this function goes here
%   Detailed explanation goes here


%% Use methods(fitresult), or methods(gof) and doc to understand what data 
% is in these results provided by curve fitting app.

%% Create Tables from fitresult and gof
fit_results_table = array2table(coeffvalues(fitresult), 'VariableNames', coeffnames(fitresult));

gof_results_table = struct2table(gof);

% Select which Columns of the table to pass back in fit_parameters_row
%fit_parameters_row = [fit_results_table, gof_results_table(:,1:2)];
fit_parameters_row = [fit_results_table, gof_results_table(:,{'sse','rsquare'})];


end

