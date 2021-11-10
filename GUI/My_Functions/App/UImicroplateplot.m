function UImicroplateplot(hAxes,assaydata)

if nargin == 1
   assaydata = hAxes; 
   hAxes = axes;
end

%% Load white to red colormap
% load microPlateAssay whiteToRed
% load whiteToBlue
cmap = parula; %whiteToBlue; %parula;
colormap(hAxes,cmap);
nbColors = size(cmap,1);

%% Color missing data in black
cmin = min(assaydata(:));
cmax = max(assaydata(:));

index = fix((assaydata(:)-cmin)/(cmax-cmin)*nbColors)+1;
index(index<1) = 1;
index(index>nbColors) = nbColors; 

rgb = ind2rgb(index,cmap);
rgb2(:,:) = rgb(:,1,:);
rgb2(isnan(assaydata(:)),:) = 0;


%% Plot wells
[nbY, nbX] = size(assaydata);
[y,x] = ind2sub([nbY,nbX], 1:nbY*nbX);

hAxes.Units = 'pixels';
sz = hAxes.Position(4);  % size=300 points is good for a axes of width=300 pixels

scatter(hAxes,x,y,sz,rgb2,'filled',...
    'MarkerEdgeColor',[0 0 0],'LineWidth',0.5);


%% Format axes
hAxes.XTick = 1:nbX;
hAxes.YTick = 1:nbY;
hAxes.TickLength = [0 0];
hAxes.YDir = 'reverse';
hAxes.XLim = [0,nbX+1];
hAxes.YLim = [0,nbY+1];
box(hAxes,'on');


%% Change YLabel to alphabetical (code from original microplateplot function)
numRows = size(assaydata,1);
if isa(assaydata,'bioma.data.DataMatrix')
    rowlabels = data.RowNames;
else
    % These are alphabetic, A-Z unless we need AA-ZZ or AAA-ZZZ,...
    
    alphaLabels = cellstr(('A':'Z')');
    if numRows <= 26
        rowlabels = alphaLabels(1:numRows);
    elseif numRows >26 && numRows <= 26^2
        rowlabels = strcat(repmat(alphaLabels,1,26)',repmat(alphaLabels,1,26));
        rowlabels = rowlabels(1:numRows);
    else
        % if we have more than 26^2 rows then we don't try to label them
        rowlabels = repmat({''},size(assaydata));
    end
end
hAxes.YTickLabel = rowlabels;


%% add colorbar
colorbar(hAxes);
caxis(hAxes,[cmin,cmax]);

