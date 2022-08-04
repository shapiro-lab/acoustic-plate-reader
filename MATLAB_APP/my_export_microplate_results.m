function [ topleads, finaltable ] = my_export_microplate_results( finaltable )
%my_export_microplate_results  exports the results to a spreadsheet
%
%   This sorts the rows of the final table by EC50 and then by row 23,
%   writes to the second tab of a a spreadsheet, then selects the top 13%
%   (approximately one top compound per plate) and writes this to the first
%   tab of a spreadsheet.
%%
finaltable = sortrows(finaltable,'EC50');
finaltable = sortrows(finaltable, -23);

writetable(finaltable,'Final_Microplate_Results.xlsx','Sheet',2);

topleads = finaltable(1:(0.13*height(finaltable)) , {'FileName', 'CompoundNumber', 'EC50', 'Reliable'});
writetable(topleads,'Final_Microplate_Results.xlsx','Sheet',1);


end

