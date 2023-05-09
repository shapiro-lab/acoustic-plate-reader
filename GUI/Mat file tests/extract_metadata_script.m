function full_map = extract_metadata(MetadataPathName, MetadataFileName)
%This function will extract information for each well from the
%metadata. It will also raise errors and return a blank if
%there are incompatability issues with the metadata file

%It will return the full_map which maps well locations to
%another map that maps parameters to values
% example full_map {A1: {color: blue, size: large, weight: 100}, A2: {color: red, size: small, weight:50}}


    function raise_error(message)
        %This is a helper function that opens up a new window and displays error messages to the user
        h=errordlg(message,'Error');
        set(h, 'WindowStyle', 'modal');
        uiwait(h);
    end

metadata_cell = readcell([MetadataPathName '\' MetadataFileName]);
size_cell = size(metadata_cell);
wells_used = metadata_cell(2:size_cell(1)); %take all the elements in the well_used column
full_map = containers.Map;

%check that each of the wells are unique
if length(wells_used) ~= length(unique(wells_used))
    raise_error('There exists a duplicate well');
    return;
end

transpose_metadata = metadata_cell'; %take transpose of metadata cell. Makes the code easier
parameters = transpose_metadata(2:size_cell(2)); %take all the parameters

%iterate through the wells, check each one, and populate full map
%The keys for full_map are strings for the wells (ex. A1), the
%values will be maps
start_index = size_cell(2)+2;
for well = wells_used
    well_string = well{1}; %Have to do this so you have a string instead of a cell containing a string
    error_msg = [well_string 'is an invalid well name'];
    
    %check if the well string starts with letters A - H
    if ~ismember(well_string(1), 'ABCDEFGH')
        raise_error(error_msg);
        full_map = containers.Map;
        return;
    end
    
    %check if the well strings ends with numbers 1 - 12
    column_num = str2double(well_string(2:end));
    if isnan(column_num)
        raise_error(error_msg);
        full_map = containers.Map;
        return
    end
    if ~ismember(column_num, 1:1:12)
        raise_error(error_msg);
        full_map = containers.Map;
        return
    end
    
    %populate the map
    well_values = transpose_metadata(start_index:start_index + size_cell(2)-2);
    full_map(well_string) = containers.Map(parameters, well_values);
end
end