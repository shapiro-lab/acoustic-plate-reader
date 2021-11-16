function [ plate ] = my_file_metadata( data, filename )
%UNTITLED11 Summary of this function goes here
%   Detailed explanation goes here

%%
% Create a table from the columns of the variable named data (a double) and add VariableNames (Column Headers)
plate = array2table(data);
plate.Properties.VariableNames = {'NegControl','Conc1','Conc2','Conc3',...
    'Conc4','Conc5','Conc6','Conc7','Conc8','Conc9','Conc10','PosControl'};
numRows = height(plate);
%%
% add filename
[~,name,ext] = fileparts(filename);
fname = string([name,ext]);
plate.FileName = repmat(fname, numRows, 1);

% add compound number
% filenum = name(end-3:end);   % or more generic (works for any number of digits at the end of file name): regexp(name,'(?<=\D+)\d+','match','once')
filenum = extractBetween(fname,'data','.csv');
filenum = double(filenum);  % convert to numeric
plate.CompoundNumber = (filenum-1)*numRows + (1:numRows)';
%%
% rearrange columns
plate = plate(:,[end-1, end, 1:end-2]);

end

