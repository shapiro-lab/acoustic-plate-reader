clear all

%% Inputs
% file parameters
pathName = 'G:\.shortcut-targets-by-id\1OPaP1L4hYKwhaXe5SqxupI4aIMb0Npgj\Acoustic plate reader project\Example data\';
SampleName = 'vramp';
saveName = pathName;

% scan_type = 'pre_post'; %'voltage_ramp', 'collapse_ramp' %TODO make these change what types of plots get made

% data parameters
disp_depth = [3 8];
noise_slice = [2 3];
disp_crange = [40 -3];
imgMode = 1; % 1 for ramping voltage, 2 for imaging voltage
computeDiff = 1; % 1 or 0 to compute pre-post-collapse difference image or not
PlateSize = [8 12]; % rows and columns of wells scanned
compar = [24 25]; % indices of voltages to compare for pre-post-collapse difference

%%
% Call raw2imgs script
raw2imgs;
load(fullfile(saveName, SampleName, 'imgs.mat'));

if imgMode == 1
    colorm = 'hot';
    FigTitle = [SampleName '_xAM-Imgs'];
elseif imgMode == 2
    colorm = 'bone';
    FigTitle = [SampleName '_Bmode-Imgs'];
end
%colorm = 'parula';

% get Z indices of sample based on disp_depth
iZd = find(Zi>=disp_depth(1)):find(Zi>=disp_depth(2));
% get Z indices of noise based on noise_slice
iZd_noise = find(Zi>=noise_slice(1)):find(Zi>=noise_slice(2));
% Convert interpolated difference image to dB, cropping based on disp_depth
Im_disp = 20*log10(squeeze(abs(dImi(iZd,:,imgMode,:))));
% Calculate noise value
noise_disp = max(mean(mean(20*log10(squeeze(abs(dImi(iZd_noise,:,imgMode,:)))))))- min(min(min(Im_disp)));
% Subtract smallest value from all values
Im_disp = Im_disp - min(Im_disp,[],[1:3]);

% get x and z sizes and number of wells
[zsize,~,Nw] = size(Im_disp);
xsize = round(zsize*6.4/diff(disp_depth));
% resize images and add RGB colormaps
% Imgs dimensions: zs, xs, colors, wells
dImgs = nan(zsize,xsize,3,Nw);
for well = 1:Nw
    % resize image
    temp = imresize(Im_disp(:,:,well),[zsize xsize]);
    % convert resized grayscale image to RGB and add to array
    dImgs(:,:,:,well) = real2rgb(temp,colorm,[noise_disp max(max(temp))+disp_crange(2)]);
end

% Plot dB difference image
prepostim = figure('Position',[877 1 561 800]);
colormap hot
montage(dImgs,'Size',PlateSize);
colorbar('Ticks',[0 1],'TickLabels',[noise_disp max(max(temp))+disp_crange(2)]);
title(FigTitle);
% savefig(fullfile(pathName, SampleName, FigTitle))
%save([saveName FigTitle],'Imgs','Im_disp');

% s = sliceViewer(Imi(:,:,:,1,1), 'Colormap',hot(256), 'ScaleFactors',[7 1 1], 'DisplayRange',[50 600]);
%% Slice and resize images and convert from grayscale to RGB
% slice Imi based on disp_depth
% Imi dimensions: zs, xs, pressures, imaging mode, wells
Imi_sliced = Imi(iZd,:,:,:,:);

% get z size and numbers of pressures, modes, and wells
[zsize,~,Np,Nm,Nw] = size(Imi_sliced);
xsize = round(zsize*6.4/diff(disp_depth));

% downsample image resolution to make size more manageable
zsize = round(zsize/3);
xsize = round(xsize/3);

% resize images and convert from grayscale to RGB
% Imgs dimensions: zs, xs, pressures, imaging modes, wells
Imgs = nan(zsize,xsize,Np,Nm,Nw);
% Im_RGB dimensions: zs, xs, colors, pressures, imaging modes, wells
Imgs_RGB = nan(zsize,xsize,3,Np,Nm,Nw);
for pressure = 1:Np
    for well = 1:Nw
        for mode = 1:Nm
            % get image from array
            im = squeeze(Imi_sliced(:,:,pressure,mode,well));
            % resize image
            im_resized = imresize(im,[zsize xsize]);
            % add resized image to array
            Imgs(:,:,pressure,mode,well) = im_resized;
            % convert grayscale image to RGB and add to array
%             Imgs_RGB(:,:,:,pressure,mode,well) = cat(3, im_resized, zeros(size(im_resized)), im_resized);
            % real2rgb changes the image scaling--do not use for
            % quantification, unmixing, etc.
            Imgs_RGB(:,:,:,pressure,mode,well) = real2rgb(im_resized,colorm,[noise_disp max(max(im_resized))+disp_crange(2)]);
        end
    end
end

% visualize one montage per pressure with sliceViewer
prepostmontages = []; % initalize montage images array
for pressure = 1:Np
    m = montage(squeeze(Imgs_RGB(:,:,1,pressure,1,:)),'Size',PlateSize); % create a montage image
    prepostmontages = cat(3,prepostmontages,m.CData); % add montage image to montage array
end
sliceViewer(prepostmontages, 'Colormap',hot(256));

% plot images across voltages
mat = cat(3,squeeze(Imgs_RGB(:,:,1,8:end-1,1,1)),squeeze(Imgs_RGB(:,:,1,8:end-1,1,3)),squeeze(Imgs_RGB(:,:,1,8:end-1,1,4)));
montage(mat,'Size',[3 14])
colormap hot
colorbar
%% Draw, quantify, and visualize ROIs
% plot image for ROI selection
imagesc(Imgs(:,:,end-1,2,1),[30 800])
colormap hot
axis image

% Draw first sample ROI
title('Click & drag to choose the sample ROI')
sample_ROI = drawrectangle;
title('Adjust ROI, then press any key to continue')
pause

% Get position of first sample ROI
sample_position = round(sample_ROI.Position);
sample_xrange = sample_position(1):sample_position(1)+sample_position(3);
sample_zrange = sample_position(2):sample_position(2)+sample_position(4);

% plot and quantify ROIs
% Quant_ROIs dimensions: pressures, imaging modes, wells
Quant_ROIs = nan(Np,2,Nw);
% Imgs_ROIs_RGB dimensions: zs, xs, colors, pressures, imaging modes, wells
Imgs_ROIs_RGB = Imgs_RGB;
for pressure = 1:Np
    for well = 1:Nw
        for mode = 1:2
            % get image from array
            im = squeeze(Imgs(:,:,pressure,mode,well));
            % create mask for ROI
            ROImask = createMask(sample_ROI);
            % quantify signal in ROI
            Quant_ROIs(pressure,mode,well) = mean(nonzeros(im.*ROImask), 'all', 'omitnan');
            % overlay blue mask for ROI on RGB image and add to array
            Imgs_ROIs_RGB(:,:,3,pressure,mode,well) = ROImask;
        end
    end
end
close all
% reshape and permute ROI quants
% Quant_ROIs new dimensions: well rows, well columns, pressures, imaging modes
Quant_ROIs = permute(reshape(Quant_ROIs, Np, 2, PlateSize(2), PlateSize(1)), [4 3 1 2]);

% view ROIs
montage(squeeze(Imgs_ROIs_RGB(:,:,:,end-1,1,:)),'Size',PlateSize);

%% Plot ROI quantifications
% plot ROI quants with microplateplot
mpplot = microplateplot(Quant_ROIs(:,:,end-1,1));
mpplot

% scan_type matters
% plot pressure vs. signal for each well in a column
legendstr = {};
% Plot ROI CNRs
fig = figure; colormap hot
voltages = [P.seed max(P.seed)+1];
columns = [1 2 3];
for column = columns
    subplot(1, length(columns), column)
    hold on
    for row = 1:8
        plot(voltages, squeeze(Quant_ROIs(row,column,:,1)))
%         plot(voltages, rescale(squeeze(Quant_ROIs(row,column,:,1)))) % normalize before plotting
        legendstr = [legendstr, {sprintf('Well %d',row)}];
    end
    title({sprintf('xAM %d',column)}),xlabel('Transducer voltage (V)'),ylabel('xAM intensity')
%     legend(legendstr)
    xlim([1 25])
    hold off
end


%%
% plot highest concentration normalized
G82L = rescale(squeeze(Quant_ROIs(1,1,:,1)));
S50C_G82L = rescale(squeeze(Quant_ROIs(1,2,:,1)));
Serratia = rescale(squeeze(Quant_ROIs(1,3,:,1)));
figure
hold on
plot(voltages, Serratia)
plot(voltages, G82L)
plot(voltages, S50C_G82L)
title('xAM'),xlabel('Transducer voltage (V)'),ylabel('Normalized xAM signal')
legend({'Serratia', 'GvpB-G82L', 'GvpB-S50C-G82L'})
xlim([10 25])
hold off

%% linear unmixing based on xAM turn-on at different voltages

% input turn-on matrix representing the fraction of each type of GV that is producing xAM signal at each applied pressure
% this can be calculated from the image data by comparing the intensities of each sample across pressures
% rows are different types of GVs; columns are different pressures
% values don't have to be scaled from 0 to 1, but all rows should have the same min and max
colmat = rescale([G82L Serratia]');

alpha = diff(colmat,1,2); % compute differential turn-on matrix "alpha" (i.e., fraction turning on at each step)

Is = nan([zsize xsize Nw Np]); % dimensions: x, y, wells, voltages
Ss = nan([zsize xsize Nw Np-1]); % dimensions: x, y, wells, voltages
Cs = nan([zsize xsize Nw size(colmat,1)]); % dimensions: x, y, wells, GV type
for well = 1:Nw
    % unmix image using differential turn-on matrix alpha
    % Im_ROIs dimensions: zs, xs, pressures, imaging modes, wells
    [Is(:,:,well,:), Ss(:,:,well,:), Cs(:,:,well,:)] = Unmix(squeeze(Imgs(:,:,:,1,well)), alpha);
end

% view difference images
montages = []; % initalize montage images array
for species = 1:Np-1
    m = montage(squeeze(Ss(:,:,:,species)),'Size',PlateSize); % create a montage image
    montages = cat(3,montages,m.CData); % add montage image to montage array
end
sliceViewer(montages, 'Colormap',hot(256),'DisplayRange',[0 100]);

% view unmixed images
montage(squeeze(Cs(:,:,:,1)),'Size',PlateSize, 'DisplayRange',[0 500]);
montage(squeeze(Cs(:,:,:,2)),'Size',PlateSize, 'DisplayRange',[0 500]);


% plot desired wells unmixed
wells = [1:5];
% scale images from 0 to 1, setting negative concentration values to 0 and
% setting very large values to 1
Cs_rescaled = rescale(Cs, 'InputMin',0, 'InputMax',max(Cs,[],[1:4])*.25);
fig = figure;
for Wi = 1:length(wells)
    subplot(3, length(wells), Wi)
    im = squeeze(Cs_rescaled(:,:,wells(Wi),1));
    % convert first image to red and plot
    imagesc(cat(3, im, zeros(size(im)), zeros(size(im))));
    set(gca,'xticklabel',[], 'yticklabel',[])
    axis square
end
for Wi = 1:length(wells)
    subplot(3, length(wells), Wi+length(wells))
    im = squeeze(Cs_rescaled(:,:,wells(Wi),2));
    % convert second image to green and plot
    imagesc(cat(3, zeros(size(im)), im, zeros(size(im))));
    set(gca,'xticklabel',[], 'yticklabel',[])
    axis square
end
for Wi = 1:length(wells)
    subplot(3, length(wells), Wi+2*length(wells))
    im1 = squeeze(Cs_rescaled(:,:,wells(Wi),1));
    % convert first image to red
    im1 = cat(3, im1, zeros(size(im1)), zeros(size(im1)));
    im2 = squeeze(Cs_rescaled(:,:,wells(Wi),2));
    % convert second image to green
    im2 = cat(3, zeros(size(im2)), im2, zeros(size(im2)));
    % sum images and plot
    imagesc(im1+im2);
    set(gca,'xticklabel',[], 'yticklabel',[])
    axis square
end

%% function definitions
function [I, D, C] = Unmix(vramp, alpha)
    % Function to unmix image voltage ramp using differential xAM turn-on matrix
    % Inputs:
    %   vramp: raw image collected at varying voltages (dimensions: x, y, voltage)
    %   alpha: differential xAM turn-on matrix
    % Outputs:
    %   I: spatially averaged raw images
    %   D: difference images
    %       for a given pixel, D contains the three observed difference signals
    %       i.e. pre minus first collapse, first collapse minus second collapse, etc
    %   C: unmixed images
    %       for a given pixel, C contains the concentration of each species
    %       C = alpha'\S; %concentrations of the three differently-collapsing species
    %       Cr = C./sum(C); %relative concentrations
    
    f = fspecial('average',[10 10]); % create averaging filter of 10x10 pixels
    I = imfilter(vramp,f); % apply spatial averaging filter to raw image
    
    % Iterating over each pixel, make difference and unmixed images
    for x = 1:size(I,1)
        for y = 1:size(I,2)
            D(x,y,:) = diff(I(x,y,:)); % get subtraction signals
            C(x,y,:) = alpha'\permute(D(x,y,:), [3 2 1]); % concentrations of the differently-collapsing species
%             C(x,y,:) = lsqnonneg(alpha',permute(D(x,y,:), [3 2 1])); % concentrations of the differently-collapsing species
        end
    end
end







%%
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