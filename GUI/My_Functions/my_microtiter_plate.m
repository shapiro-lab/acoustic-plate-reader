function [ plate_result ] = my_microtiter_plate( FileName, filenum )
%MY_MICROTITER_PLATE completes entire analysis of a single 96-well plate

% This function consists of several custom functions:
% 1. Import data from .csv files of microwell plate.  The example is 
% set to 96 well data (8 rows by 12 columns) but is user-configurable 
% inside the my_import_plate function.

% 2. Perform positive and negative controls on each file

% 3. Curve fit to a four parameter model

% 4. Quality control on the curve fit data

%% Custom function to load data, add Variable Headers, add File Name and hand back a table
plate = my_import_plate(FileName{filenum});                       

%% Custom function to apply positive and negative controls to the table and overwrite new results to the old table
plate = my_plate_controls(plate);                                       

% iterate through each row of one plate: add curve fitting and qc data
for ii = 1:height(plate)                                                
    % assign to Y the preprocessed data (change X inside the function if needed)
    Y = plate{ii,4:13};                                                 
    
    % CURVE FITTING
    % Custom function: pass in Y to our custom 4 parameter fit, get out a table row
    fitting = my_four_parameter_fit(Y);                                 
    % build up the platefit table table of eight rows (one plate)
    plate_fit(ii,:) = fitting;                                          
    % $$ y = min + (max-min)/(1 + (10^x/EC50)^(Hills_Slope)) = 0$$
    
    
    % QUALITY CONTROL
    % Custom function: pass EC50 and Hills_Slope into the QC function and get back a QC table row
    qc = my_qc(fitting.EC50, fitting.Hills_Slope);                      
    % build up the Quality Control (qc) table
    plate_qc(ii,:) = qc;                                                
    
end

% build the result for each plate to be passed back out of the function
plate_result = [plate, plate_fit, plate_qc];

end

