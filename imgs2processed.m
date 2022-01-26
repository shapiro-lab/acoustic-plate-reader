clear all

%% Inputs
% file parameters
pathName = 'C:\Users\Administrator\Dropbox\GV Team\verasonics system\Vantage-4.6.2-RCH\Data\';
SampleName = '220126';

% scan_type = 'pre_post'; %'voltage_ramp', 'collapse_ramp' % TODO make these change what types of plots get made

% data parameters
disp_crange = [40 -3]; % limits of colorbar
imgMode = 1; % 1 for ramping voltage, 2 for imaging voltage
computeDiff = 1; % 1 or 0 to compute pre-post-collapse difference image or not
compar = [1 2]; % indices of voltages to compare for pre-post-collapse difference
trans = 'L22'; % L22 or L10

%%
% specify transducer parameters
if trans == 'L22'
    disp_depth = [3 10];
    noise_slice = [2 3];
    ScaleFactors = [1 1 1];
elseif trans == 'L10'
    disp_depth = [18 25];
    noise_slice = [10 12];
    ScaleFactors = [3 1 1];
end

% Call raw2imgs script
raw2imgs;
load(fullfile(pathName, SampleName, 'imgs.mat'));

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
savefig(fullfile(pathName, SampleName, FigTitle))

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
rawMontages = []; % initalize montage images array
for pressure = 1:Np
    m = montage(squeeze(Imgs(:,:,pressure,1,:)),'Size',PlateSize); % create a montage image
    rawMontages = cat(3,rawMontages,m.CData); % add montage image to montage array
end
sliceViewer(rawMontages, 'Colormap',hot(256), 'ScaleFactors',ScaleFactors, 'DisplayRange',[20 400]);

% visualize one montage per pressure with sliceViewer
scaledMontages = []; % initalize montage images array
for pressure = 1:Np
    m = montage(squeeze(Imgs_RGB(:,:,1,pressure,1,:)),'Size',PlateSize); % create a montage image
    scaledMontages = cat(3,scaledMontages,m.CData); % add montage image to montage array
end
sliceViewer(scaledMontages, 'ScaleFactors',ScaleFactors, 'Colormap',hot(256));

% plot images across voltages
% Im_RGB dimensions: zs, xs, colors, pressures, imaging modes, wells
mat = cat(3,squeeze(Imgs_RGB(:,:,1,:,1,1)),squeeze(Imgs_RGB(:,:,1,:,1,2)),squeeze(Imgs_RGB(:,:,1,:,1,3)));
montage(mat,'Size',[3 Np])
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
colormap hot
colorbar
mpplot
savefig(fullfile(pathName, SampleName, [FigTitle '_ROIs']))

% scan_type matters
% plot pressure vs. signal for each well in a column
legendstr = {};
% Plot ROI CNRs
figure; colormap hot
voltages = [P.seed max(P.seed)+1];
columns = [1:4];
for column = columns
    subplot(1, length(columns), find(columns==column))
    hold on
    for row = 1:2
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
% plot all samples of a given type
% Quant_ROIs dimensions: well rows, well columns, pressures, imaging modes
% pressure_range = voltages;
pressure_range = [5:length(voltages)-1];
Serratia = reshape(squeeze(Quant_ROIs(1:2,1,pressure_range,1)), [], length(pressure_range));
Q65Sfs = reshape(squeeze(Quant_ROIs(1:2,2,pressure_range,1)), [], length(pressure_range));
A68R = reshape(squeeze(Quant_ROIs(1:2,3,pressure_range,1)), [], length(pressure_range));
B230 = reshape(squeeze(Quant_ROIs(1:2,4,pressure_range,1)), [], length(pressure_range));
% Mega13 = reshape(squeeze(Quant_ROIs(4,:,:,1)), [], length(pressure_range));
% Mega14 = reshape(squeeze(Quant_ROIs(5,:,:,1)), [], length(pressure_range));
% Serratia = reshape(squeeze(Quant_ROIs(6,:,:,1)), [], length(pressure_range));


% % plot all replicates of each sample
figure;
hold on
plot(voltages(pressure_range), reshape(A68R,[],length(pressure_range)))
plot(voltages(pressure_range), reshape(Q65Sfs,[],length(pressure_range)))
plot(voltages(pressure_range), reshape(Serratia,[],length(pressure_range)))
plot(voltages(pressure_range), reshape(B230,[],length(pressure_range)))
% plot(voltages, reshape(AC,[],length(voltages)))
% plot(voltages, reshape(Serratia,[],length(voltages)))
% plot(voltages, reshape(S50C_G82L,[],length(voltages)))
title('xAM'),xlabel('Transducer voltage (V)'),ylabel('xAM signal')
legend({'A68R', 'Q65Sfs', 'Serratia', 'GvpB-S50C-G82L'})
% xlim([2 25])
hold off

% plot average of each sample
figure;
hold on
plot(voltages(pressure_range), mean(A68R))
plot(voltages(pressure_range), mean(Q65Sfs))
plot(voltages(pressure_range), mean(Serratia))
plot(voltages(pressure_range), mean(B230))
% plot(voltages(pressure_range), mean(Mega13))
% plot(voltages(pressure_range), mean(Mega14))
% plot(voltages(pressure_range), mean(Serratia))
title('xAM'),xlabel('Transducer voltage (V)'),ylabel('xAM signal')
legend({'A68R', 'Q65Sfs', 'Serratia', 'GvpB-S50C-G82L'}, 'Location','northwest', 'Interpreter', 'none')
% legend({'pMetTU1-A_Sv6K_lacI-Ptac-lacO_RC_AnaACNJKFGVW_BBa-B0015','pMetTU1-A_Sv7K_araC-PBAD_RC_AnaA_T_AnaCNJKFGW','pMetTU1-A_Sv5K_lacI-Ptac-lacO_RC_AnaA-A68R_T_AnaC-MegaRNFGLSKJTU_Bba-B0015','Mega13','Mega14','Serratia'}, 'Location','northwest', 'Interpreter', 'none')
% xlim([10 20])
set(findall(gca, 'Type', 'Line'),'LineWidth',2);
hold off

figure;
hold on
plot(voltages(pressure_range), mean(Ana1))
plot(voltages(pressure_range), mean(Ana2))
plot(voltages(pressure_range), mean(AC))
plot(voltages(pressure_range), mean(Mega13))
title('xAM'),xlabel('Transducer voltage (V)'),ylabel('xAM signal')
legend({'Sv6K_Ptac-lacO_AnaACNJKFGVW_BBa-B0015','Sv7K_PBAD_AnaACNJKFGW','Sv5K_Ptac-lacO_AnaA-A68R_AnaC-MegaRNFGLSKJTU_Bba-B0015','Mega13'}, 'Location','northwest', 'Interpreter', 'none')
xlim([10 20])
set(findall(gca, 'Type', 'Line'),'LineWidth',2);
hold off

%%
% plot highest concentration normalized
% Quant_ROIs dimensions: well rows, well columns, pressures, imaging modes
% Serratia = rescale(mean(squeeze(Quant_ROIs(1:2,1,pressure_range,1))));
Q65Sfs = rescale(mean(squeeze(Quant_ROIs(1:2,2,pressure_range,1))));
A68R = rescale(mean(squeeze(Quant_ROIs(1:2,3,pressure_range,1))));
B230 = rescale(mean(squeeze(Quant_ROIs(1:2,4,pressure_range,1))));

figure;
hold on
% plot(voltages(pressure_range), Serratia)
plot(voltages(pressure_range), Q65Sfs)
plot(voltages(pressure_range), A68R)
plot(voltages(pressure_range), B230)
title('xAM (Normalized)'),xlabel('Transducer voltage (V)'),ylabel('Normalized xAM signal')
legend({'Q65Sfs', 'A68R', 'GvpB-S50C-G82L'})
% legend({'Sv5K_Ptac-lacO_AnaA-A68R_AnaC-MegaRNFGLSKJTU_Bba-B0015', 'Serratia', 'GvpB-S50C-G82L'}, 'Location','northwest', 'Interpreter', 'none')
% xlim([5 20])
set(findall(gca, 'Type', 'Line'),'LineWidth',2);
hold off

%% linear unmixing based on xAM turn-on at different voltages

% input turn-on matrix representing the fraction of each type of GV that is producing xAM signal at each applied pressure
% this can be calculated from the image data by comparing the intensities of each sample across pressures
% rows are different types of GVs; columns are different pressures
% values don't have to be scaled from 0 to 1, but all rows should have the same min and max
colmat = rescale([A68R' B230']');

alpha = diff(colmat,1,2); % compute differential turn-on matrix "alpha" (i.e., fraction turning on at each step)

Is = nan([zsize xsize Nw length(pressure_range)]); % dimensions: x, y, wells, voltages
Ss = nan([zsize xsize Nw length(pressure_range)-1]); % dimensions: x, y, wells, voltages
Cs = nan([zsize xsize Nw size(colmat,1)]); % dimensions: x, y, wells, GV type
for well = 1:Nw
    % unmix image using differential turn-on matrix alpha
    % Imgs dimensions: zs, xs, pressures, imaging modes, wells
    [Is(:,:,well,:), Ss(:,:,well,:), Cs(:,:,well,:)] = Unmix(squeeze(Imgs(:,:,pressure_range,1,well)), alpha);
end

% view sequential difference images
montages = []; % initalize montage images array
for pressure = 1:length(pressure_range)-1
    m = montage(squeeze(Ss(:,:,:,pressure)),'Size',PlateSize); % create a montage image
    montages = cat(3,montages,m.CData); % add montage image to montage array
end
figure;
sliceViewer(montages, 'Colormap',hot(256),'DisplayRange',[0 100]);
title('Sequential difference images')

% view unmixed images
figure;
montage(squeeze(Cs(:,:,:,1)),'Size',PlateSize, 'DisplayRange',[0 500]);
title('Species 1')
figure;
montage(squeeze(Cs(:,:,:,2)),'Size',PlateSize, 'DisplayRange',[0 500]);
title('Species 2')


% plot desired wells (indexed across rows) unmixed
wells = [5 8];
% scale images from 0 to 1, setting negative concentration values to 0 and
% setting very large values to 1
Cs_rescaled = rescale(Cs, 'InputMin',0, 'InputMax',max(Cs,[],[1:4])*.25);
figure;
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
    im = squeeze(rescale(Cs_rescaled(:,:,wells(Wi),2)));
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

%%

Vs = nan([Np-1 Np-1 Nw]); % dimensions: voltages, voltages, wells
Ss = nan([Np-1 Np-1 Nw]); % dimensions: voltages, voltages, wells
USs = nan([zsize xsize Np-1 Nw]); % dimensions: x, y, voltages, wells

for well = 1:1
    % run SVD on each well
    % Im_ROIs dimensions: zs, xs, pressures, imaging modes, wells
    [Vs(:,:,well), Ss(:,:,well), USs(:,:,:,well)] = svdNplot(squeeze(Imgs(:,:,1:end-1,1,well)));
end

%% function definitions
function [I, D, C] = Unmix(vramp, alpha)
    % Function to unmix image voltage ramp using differential xAM turn-on matrix
    % Inputs:
    %   vramp: raw image collected at varying voltages (dimensions: x, y, voltage)
    %   alpha: differential xAM turn-on matrix
    % Outputs:
    %   I: spatially averaged raw images
    %   D: sequential difference images
    %       for a given pixel, D contains the calculated sequential difference signals
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
%             size(D)
%             size(permute(D(x,y,:), [3 2 1]))
%             size(alpha')
            C(x,y,:) = alpha'\permute(D(x,y,:), [3 2 1]); % concentrations of the differently-collapsing species
%             C(x,y,:) = lsqnonneg(alpha',permute(D(x,y,:), [3 2 1])); % concentrations of the differently-collapsing species
        end
    end
end


function [V, S, US_resh] = svdNplot(Mat4svd)
    % Mat4svd is a 3D matrix with dimensions: X, Y, Voltage

    Mat_reshaped = reshape(Mat4svd,size(Mat4svd,1)*size(Mat4svd,2),size(Mat4svd,3));
    
    covar_mat = Mat_reshaped' * Mat_reshaped; % To save calculation time, you can calculate the SVD on the covariance matrix of your Mat_reshape
    
    % SVD: S = lambda x U x V
    % with lamba = "weight" of your modes
    %      U = spatial modes
    %      V = voltage modes
    [V,S] = svd(covar_mat);
    US = Mat_reshaped * V;
    US_resh = reshape(US,size(Mat4svd,1),size(Mat4svd,2),size(US,2));
    
    % Display the first 8 spatial modes
    figure; 
    for sMode = 1:8
        subplot(2,4,sMode)
        tmp = US_resh(:,:,sMode);
        imagesc(tmp);
        colormap jet
        caxis([-max(abs(tmp(:))) max(abs(tmp(:)))]);
        colorbar
        title(sMode)
    end
    
    % Display the first 4 voltage modes
    figure;
    for vMode = 1:4
        subplot(2,2,vMode)
        plot(V(:,vMode))
        title(vMode)
    end
end