overlay = 0;
imgMode = 1;
saveName = '/Users/Sanyo 1/Documents/MATLAB/Vantage-3.3.0-1710061400/Data/PlateReader/';
SizeThreshold = 300;
test = nan(1,96);
ROISizes = test;
h = waitbar(0,'Processing...');
for i = 1:total_n
Ind = i;
autoROI1;
test(i) = autoROI;
ROISizes(i) = ROISize;
waitbar(i/total_n);
if overlay
    pause(0.5);
    close(fig);
end
end
close(h);
test1 = reshape(test,12,8);
ROISize1 = reshape(ROISizes,12,8);
test1 = test1';
ROISize1 = ROISize1';
name1 = reshape(PlateCoordinate,12,8);
name1 = name1';
name2 = reshape(1:96,12,8);
name2 = name2';
test2 = (test1 ./ ROISize1) .* (ROISize1 > SizeThreshold);
if imgMode == 1
    FigTitle = [SampleName 'xAM'];
elseif imgMode == 2
    FigTitle = [SampleName 'Bmode'];
end
fig = PlatePlot1(test2,FigTitle);
savefig([saveName FigTitle])
save([saveName FigTitle],'test1','ROISize1','test2','name1','name2');
function fig_out = PlatePlot1(data,fig_title,fontsize,varargin)
if nargin == 2
    fontsize = 16;
end
fig_out = figure;
[y,x] = size(data);
imagesc(data);axis image
xticks(1:x)
yticks(1:y)
ylabelnames = {'A','B','C','D','E','F','G','H'};
ylabelnames = ylabelnames(1:y);
yticklabels(ylabelnames)
colorbar;
title(fig_title)
set(gca,'fontsize',fontsize);
end