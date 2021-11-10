function data = my_import_file(filename)
display(class(filename));
%IMPORTFILE Import numeric data from a text file as a matrix.
%   MICROTITERDATA0001 = IMPORTFILE(FILENAME) Reads data from text file
%   FILENAME for the default selection.
%
%% Post processing for unimportable data.
% No unimportable data rules were applied during the import, so no post
% processing code is included. To generate code which works for
% unimportable data, select unimportable cells in a file and regenerate the
% script.

%% Create output variable
data = readmatrix(filename);