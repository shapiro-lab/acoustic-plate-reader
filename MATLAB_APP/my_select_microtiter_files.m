function [FileNames, PathName] = my_select_microtiter_files() 
%Select which files to use
%   To prototype we used autogenerated code to load one file from a
%   directory.  

%   This function uses uigetfile to allow the user to tell MATLAB a number 
%   of files to use at once.

%% What files to load?

% Point MATLAB to your *.xlsx files in a directory 
[FileNames,PathName] = uigetfile({'*.csv';'*.xlsx';'*.xls'},'Select Data Files to Use','MultiSelect','on'); 

% if only one file was selected, convert the class of FileName from char to cell 
if ischar(FileNames) 
    FileNames={FileNames}; 
end          

% % we could use this later in the app if we want to display the FullFileName
% FullFileName = fullfile(PathName,FileNames)

end