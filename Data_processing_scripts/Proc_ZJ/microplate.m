% Load a MAT-file, included with the Bioinformatics Toolbox? software, which contains two variables: assaydata, an 8-by-12 matrix of data values from a microtiter plate, and whiteToRed, a 64-by-3 matrix that defines a colormap.
load microPlateAssay

% % Overlay an X on well E8.
% Create an empty cell array.
mask = cell(8,12);
% Add the string 'X' to the cell in the fifth row and eighth column of the array.
mask{5,8} = 'X';
% Pass the cell array to the microplateplot function using the 'TextLabels' property.
microplateplot(assaydata,'TEXTLABELS',mask);
    
% Add colorbar to plot
colorbar

% Change the visualization to use a hot colormap, and then view a tooltip displaying the value of any well by clicking and holding on the well.
colormap(hot)
% Notice that all wells in column 12 are black, indicating missing data.