function [ pre_processed_data ] = my_plate_controls(data)
% PREPROCESS OUR DATA
%
% We need to do two operations to our raw data - a negative control and the
% positive control - and pass back a pre-processed array for further
% analysis
%
%% Apply the negative and positive controls
% NEGATIVE CONTROL: Subtract the mean of the first column of experimental data from the value in each well in all the wells
data = data - mean(data(:,1));

% POSITIVE CONTROL: Scale (normalize) each row so that the Postitive Control (in column 12, or the last value) is 100
data = data./data(:,end)*100;      

% overwrite the original values of the raw experimental data in the table called plate
pre_processed_data = data;
    
end

