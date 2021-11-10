function [ plate_result ] = my_A_Microtiter_Plate_Onefile( filename )

%MY_ONE_PLATE This is literally cut and pasted from the file
% A_Microtiter_Plate_Onefile, except we turned that script into a fucnction
% by adding a function call with arguments at the top row and and end in
% the last row of this file.

%% Microtiter Plate Analyisis - Single File Workflow
% This file completes the entire analysis of a single 96-well plate file
%
% 1. Import one.csv file from an instrument into a varable (a double)
%
% 2. Apply positive and negative controls on each file
% 
% 3. Add file and compound metadata to a table
%
% 4. Curve fit to a four parameter model to sigmoidal curve and fitresults
%
% 5. Turn fitresults into a table for each row
%
% 6. Quality control on the curve fit data from fitresults
%
% 7. Merge the subtables generated at each step into a single table
%
%% LOAD DATA - comment out for the app or UIselect
%filename = 'microtiter_data0001.csv';

% A Custom function to load data

% check data type and convert
if ~ischar(filename)
    filename = char(filename);
end

data = my_import_file(filename);

%% VISUALIZE THE RAW DATA (Comment out when you run a lot of files)
% figure
% bar3(data)
% 
% figure
% boxplot(data)
% 
% figure
% waterfall(data)
% 
% figure
% microplateplot(data)
%% APPLY POSITIVE AND NEGATIVE CONTROLS
data = my_plate_controls(data);                                      

%% ASSOCIATE FILE NAMES AND COMPOUND NAMES
plate = my_file_metadata(data, filename);

%% ITERATE THROUGH EACH ROW OF THE PLATE FOR CURVE FIT AND QC DATA
% let's say 1 <---> 1 uM so we're ranging from 1 nM to 1 M, for example
X = log10([0.001, 0.01, 0.1, 1, 10, 100, 1000, 10000, 100000, 1000000]);

for ii = 1:height(plate)                                                
    Y = plate{ii,4:13};  
    
    % CURVE FIT: pass a preprocessed row of data (Y) to our custom 4
    % parameter fit and get out a curve and associated parameters
    [fitresult, gof] = my_four_parameter_fit(X, Y );                               
    % build up the platefit table table of eight rows (one plate)
    [ fit_parameters_row(ii,:) ] = my_fit_parameters( fitresult, gof );
 
    % QUALITY CONTROL: pass EC50 and Hills_Slope into the this custom QC
    % function inspired by the literature and get back a QC table row
    qc = my_qc(fit_parameters_row.EC50(ii), fit_parameters_row.Hills_Slope(ii));                      
    % build up the Quality Control (qc) table
    plate_qc(ii,:) = qc;                                                
    
end

%% MERGE RESULTS FROM EACH SECTION INTO A SINGLE TABLE
plate_result = [plate, fit_parameters_row, plate_qc];

end

