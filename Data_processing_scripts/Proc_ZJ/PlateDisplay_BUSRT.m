%clear all

saveName = 'C:\Users\verasonics\Dropbox\GV Team\verasonics system\Vantage-4.6.2-RCH\Data\Wheesoo_221119_BURST/';
pathName = 'C:\Users\verasonics\Dropbox\GV Team\verasonics system\Vantage-4.6.2-RCH\Data\Wheesoo_221119_BURST/';
%ExperimentDate = '210329_RCH';
SampleName = 'GV3_2022-11-19@14-41-402';

%PlateProc1;

imgMode = 2;
disp_depth = [5 10];
noise_slice = [2 3];
disp_crange = [0 0];
if imgMode == 1
    colorm = 'hot';
else
    colorm = 'bone';
end
%colorm = 'jet';
PlateSize = [6 6];

if imgMode == 1
    FigTitle = [SampleName '-BURST-Imgs'];
elseif imgMode == 2
    FigTitle = [SampleName '-(Pre-Post)-Imgs'];
end

iZd = find(Zi>=disp_depth(1)):find(Zi>=disp_depth(2));
iZd_noise = find(Zi>=noise_slice(1)):find(Zi>=noise_slice(2));
Im_disp = squeeze(abs(dImi(iZd,:,imgMode,:)));
noise_disp = max(mean(mean((squeeze(abs(dImi(iZd_noise,:,imgMode,:))))))) - min(min(min(Im_disp)));
Im_disp = Im_disp - min(min(min(Im_disp)));
[zs,~,Ns] = size(Im_disp);
xs = round(zs/diff(disp_depth)*(Xi(end)-Xi(1)));
Imgs = nan(zs,xs,3,Ns);
for dispi = 1:Ns
    temp = imresize(Im_disp(:,:,dispi),[zs xs]);
    Imgs(:,:,:,dispi) = real2rgb(temp,colorm,[noise_disp-3 max(max(temp))+disp_crange(2)]);
end
fig = figure('Position',[877 1 561 800]);
montage(Imgs,'Size',PlateSize);
title(FigTitle);

PlatePlot1(BURST_sum,FigTitle);
%savefig([saveName FigTitle])
%save([saveName FigTitle],'Imgs','Im_disp');

% saveName = '/Users/Sanyo 1/Documents/MATLAB/Vantage-3.3.0-1710061400/Data/PlateReader/';
% SizeThreshold = 300;
% test = nan(1,96);
% ROISizes = test;
% h = waitbar(0,'Processing...');
% for i = 1:total_n
% Ind = i;
% autoROI1;
% test(i) = autoROI;
% ROISizes(i) = ROISize;
% waitbar(i/total_n);
% if overlay
%     pause(0.5);
%     close(fig);
% end
% end
% close(h);
% test1 = reshape(test,12,8);
% ROISize1 = reshape(ROISizes,12,8);
% test1 = test1';
% ROISize1 = ROISize1';
% name1 = reshape(PlateCoordinate,12,8);
% name1 = name1';
% name2 = reshape(1:96,12,8);
% name2 = name2';
% test2 = (test1 ./ ROISize1) .* (ROISize1 > SizeThreshold);
% if imgMode == 1
%     FigTitle = [SampleName 'xAM'];
% elseif imgMode == 2
%     FigTitle = [SampleName 'Bmode'];
% end
% fig = PlatePlot1(test2,FigTitle);
% savefig([saveName FigTitle])
% save([saveName FigTitle],'test1','ROISize1','test2','name1','name2');
% function fig_out = PlatePlot1(data,fig_title,fontsize,varargin)
% if nargin == 2
%     fontsize = 16;
% end
% fig_out = figure;
% [y,x] = size(data);
% imagesc(data);axis image
% xticks(1:x)
% yticks(1:y)
% ylabelnames = {'A','B','C','D','E','F','G','H'};
% ylabelnames = ylabelnames(1:y);
% yticklabels(ylabelnames)
% colorbar;
% title(fig_title)
% set(gca,'fontsize',fontsize);
% end