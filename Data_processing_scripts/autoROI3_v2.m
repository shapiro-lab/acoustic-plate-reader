close all
clear iZt ixZtemp

%%
% sample_depth = [1 9.1]; % display/sample depth range in mm
sample_depth = [1 8]; % display/sample depth range in mm
DisplayMode = 1; % How to display images: 0 = show nothing, 1 = always show images + ROIs, 2 = only show images + ROIs with low-confidence ROI selection
testfilter = 0; % 0 or 1 to display the thresholding filter
savedata = 1; % 0 or 1 to save processed data
xROI_size = 10; % size of the ROI in x dimension    
zROI_size = 40; % size of the ROI in z dimension 
filter_size1 = [40 5]; % size of the median filter in [z,x] dimension
noise_stdratio = 2; % ratio of the noise STD added to the noise mean to decide noise floor 
noise_slices = repmat([2 3],total_n,1); % default noise depth (in mm) for noise ROI selection
noise_slices_backup = [8 9;1 2;7 8]; % backup noise depth (in mm) for noise ROI selection in case of abnormal noise level
noiseThresh = 100; % threshold of noise level (mV) for abnormal high noise level / changing noise depth slices
TemplateMode = 3; % ROI template to use: 1 = filled wells, 2/3 = empty wells with thinner (2) or thicker (3) well wall, use 1 as default and 2/3 if well is not full
if ~exist('WellTemplates','var') % load the ROI template if needed, please put the .mat file in accessible folders
    load('WellTemplates.mat');
end
WellTemplate = WellTemplates{TemplateMode}; % Apply ROI template as selected.
CST = ConfidenceScoreThreshold(TemplateMode); % Apply corresponding confidence score threshold for reminder to adjust ROIs (abitrary defined thresholds)

figs = gobjects(1,total_n); % initialize array to hold figure objects
noiseZ = nan(2,total_n); % initialize array for noise slice bounds
confScore = nan(1,total_n); % initialize array for ROI prediction confidence score 
ROI_Centers = nan(total_n,2); % initialize array for locations of ROI centers
sampROI_means = nan(Nf,2,total_n); % initialize array for sample means
noiseROI_means = nan(Nf,2,total_n); % initialize array for noise means
noiseROI_stds = nan(Nf,2,total_n); % initialize array for noise ROI STDs

%% manual correction of ROI selection
xcorrection = zeros(1,total_n); % initialize array for ROI correction in x dimension (correction in the unit of voxel count, ~0.1 mm/count)
zcorrection = zeros(1,total_n); % initialize array for ROI correction in z dimension (correction in the unit of voxel count, ~0.0124 mm/count)  
skip = zeros(1,total_n); % ROI correction in well indices to skip processing

% xcorrection([2]) = 7;
% xcorrection([29]) = -5;
% xcorrection([61]) = 25;
% xcorrection([74]) = 25;
% xcorrection([76]) = 55;
% xcorrection([85]) = 25;
% xcorrection([86]) = 30;
% xcorrection([18]) = 5;
% xcorrection([34]) = -5;
% xcorrection([38]) = 5;
% xcorrection([60]) = -5;
% xcorrection([27]) = 5;
% xcorrection([53]) = -5;
% xcorrection([59]) = -5;
% xcorrection([85]) = -5;
%xcorrection([37]) =55;

% zcorrection(1:total_n) = 20;
% zcorrection(1) = -20;
% zcorrection(73) = 20;
% zcorrection(74) = 20;
% zcorrection(79) = 20;
% zcorrection(84) = 20;
% zcorrection(86) = -80;
% zcorrection(25) = 40;
% zcorrection(31) = 20;
% zcorrection(81) = -30;
% zcorrection(96) = -10;
zcorrection = zcorrection + 20;
%% quantify ROIs
for wellIx = 1:total_n
    if ~skip(wellIx)
        ImTemp = Imi{1,2,wellIx}; % Fetch Bmode image at first voltage to use for ROI detection
        noise_slice = noise_slices(wellIx,:); % Fetch noise depth slice to select noise ROI
        if ~exist('ixZtemp','var') % Select display depth if not done yet
            ixZtemp = find(Zi>=sample_depth(1),1,'first'):find(Zi>=sample_depth(2),1,'first'); 
        end
        ImTemp = ImTemp(ixZtemp,:); % Crop the images to defined sample depth to decrease data size
        ZixTemp = Zi(ixZtemp); % Define new z-axis after cropping
        noiseZ(:,wellIx) = [find(ZixTemp>=noise_slice(1),1,'first') find(ZixTemp>=noise_slice(2),1,'first')]; % Find start/end indices for corresponding noise ROI
        iZd_noise = noiseZ(1,wellIx):noiseZ(2,wellIx); % Make vector of noise Z indices
        noiseROItemp = mean(mean(ImTemp(iZd_noise,:))) + noise_stdratio * std2(ImTemp(iZd_noise,:)); % Calculate the noise level for thresholding 
        i_noiseThresh = 1; % Initialize counting for calculating based on backup noise depth slices 
        while noiseROItemp > noiseThresh && (i_noiseThresh < length(noise_slices_backup(:,1))) % Change noise depth slices until below the threshold or run out of slices
            noise_slice = noise_slices_backup(i_noiseThresh,:); % Fetch new slice
            noiseZ(:,wellIx) = [find(ZixTemp>=noise_slice(1),1,'first') find(ZixTemp>=noise_slice(2),1,'first')]; % Make array of noise Z bounds
            iZd_noise = noiseZ(1,wellIx):noiseZ(2,wellIx); % Update vector of noise Z indices
            noiseROItemp = min(noiseROItemp,mean(mean(ImTemp(iZd_noise,:))) + noise_stdratio * std2(ImTemp(iZd_noise,:))); % Take lower one between the old and new noise level 
            i_noiseThresh = i_noiseThresh + 1; % Iterate 
        end
        
        ImTemp_filt1 = medfilt2(ImTemp,filter_size1); % Apply median filter to remove salt&pepper noise
        ImTemp_filt2 = (ImTemp_filt1 > noiseROItemp); % Apply thresholding based on the noise level
        if testfilter % Debugging mode, display applied filters
            figure;
            imagesc(Xi, ZixTemp, ImTemp_filt2); axis image;
        end
        cor = xcorr2(WellTemplate,double(ImTemp_filt2)); % Calculate 2D cross-correlation between images and centered ROI template 
        [zpeak,xpeak] = find(cor==max(cor(:)),1); % Find the max correlation location
        confScore(wellIx) = cor(zpeak,xpeak); % Use the correlation at max location as confidence score
        zoffset = zpeak - size(ImTemp_filt2,1) + zcorrection(wellIx); % Convert max correlation location to index offset in z dimension and apply manual z correction
        xoffset = xpeak - size(ImTemp_filt2,2) - xcorrection(wellIx); % Convert max correlation location to index offset in x dimension and apply manual x correction
        ROI_Centers(wellIx,:) = TemplateCenter - [zoffset xoffset]; % Apply the offset to define the ROI center
        ZixROI = [ROI_Centers(wellIx,1) - zROI_size ROI_Centers(wellIx,1) + zROI_size]; % Extend from center to get start/end indices in z dimension 
        XixROI = [ROI_Centers(wellIx,2) - xROI_size ROI_Centers(wellIx,2) + xROI_size]; % Extend from center to get start/end indices in x dimension
        if ZixROI(1) < 0
            ZixROI(1) = 1;
            ZixROI(2) = 1+2*zROI_size;
            disp(['Warning: well#' num2str(wellIx) ' might be off.'])
        elseif ZixROI(2) > length(ZixTemp)
            ZixROI(2) = length(ZixTemp);
            ZixROI(1) = ZixROI(2)-2*zROI_size;
            disp(['Warning: well#' num2str(wellIx) ' might be off.'])
        end
        if XixROI(1) < 0
            XixROI(1) = 1;
            XixROI(2) = 1+2*xROI_size;
            disp(['Warning: well#' num2str(wellIx) ' might be off.'])
        elseif XixROI(2) > length(Xi)
            XixROI(2) = length(Xi);
            XixROI(1) = XixROI(2)-2*xROI_size;
            disp(['Warning: well#' num2str(wellIx) ' might be off.'])
        end
        if DisplayMode == 1 % Always-on display mode, show images and predicted ROIs
            figure(wellIx);
            imagesc(Xi, ZixTemp, 20*log10(abs(ImTemp)), [20 80]); axis image; colormap hot;
            hold on;
            noiseROI = drawrectangle('Position',[Xi(1) ZixTemp(iZd_noise(1)) Xi(end)-Xi(1) ZixTemp(iZd_noise(end))-ZixTemp(iZd_noise(1))], 'EdgeColor','w', 'LineWidth',2, 'Tag',['noise_' num2str(wellIx)]); % Draw noise ROI
            addlistener(noiseROI,'ROIMoved',@updateNoiseROI);
            sampROI = drawrectangle('Position',[Xi(XixROI(1)) ZixTemp(ZixROI(1)) Xi(XixROI(2))-Xi(XixROI(1)) ZixTemp(ZixROI(2))-ZixTemp(ZixROI(1))], 'EdgeColor','g', 'LineWidth',2, 'Tag',['samp_' num2str(wellIx)]); % Draw sample ROI
            addlistener(sampROI,'ROIMoved',@updateSampROI);
            title(['Frame #' num2str(wellIx)]);
            figs(wellIx) = figure(wellIx); % store figure object in array
            pause(0.3);
        elseif DisplayMode == 2 && (confScore(wellIx) < CST) % Only display images below confidence score threshold
            figure(wellIx);
            imagesc(Xi, ZixTemp, 20*log10(abs(ImTemp)), [20 80]); axis image; colormap hot;
            hold on;
            noiseROI = drawrectangle('Position',[Xi(1) ZixTemp(iZd_noise(1)) Xi(end)-Xi(1) ZixTemp(iZd_noise(end))-ZixTemp(iZd_noise(1))], 'EdgeColor','w', 'LineWidth',2, 'Tag',['noise_' num2str(wellIx)]); % Draw noise ROI
            addlistener(noiseROI,'ROIMoved',@updateNoiseROI);
            sampROI = drawrectangle('Position',[Xi(XixROI(1)) ZixTemp(ZixROI(1)) Xi(XixROI(2))-Xi(XixROI(1)) ZixTemp(ZixROI(2))-ZixTemp(ZixROI(1))], 'EdgeColor','g', 'LineWidth',2, 'Tag',['samp_' num2str(wellIx)]); % Draw sample ROI
            addlistener(sampROI,'ROIMoved',@updateSampROI);
            title(['Frame #' num2str(wellIx)]);
            figs(wellIx) = figure(wellIx); % store figure object in array
            pause(0.3);
        end

        for frame = 1:Nf % Calulate sample mean, noise mean, and noise STD using the predicted ROI for all frames
            for imMode = 1:2
                ImTemp = Imi{frame,imMode,wellIx}(ixZtemp,:);
                sampROI_means(frame,imMode,wellIx) = mean(mean(ImTemp(ZixROI(1):ZixROI(2),XixROI(1):XixROI(2))));
                noiseROI_means(frame,imMode,wellIx) = mean(mean(ImTemp(iZd_noise,:)));
                noiseROI_stds(frame,imMode,wellIx) = std2(ImTemp(iZd_noise,:));
            end
        end
    end
end
sampCNR = 20 * log10(abs(sampROI_means - noiseROI_means) ./ noiseROI_stds); % Calculate sample CNR
AM_Bmode_ratio = 20*log10(abs((sampROI_means(:,1,:) - noiseROI_means(:,1,:)) ./ (sampROI_means(:,2,:) - noiseROI_means(:,1,:)))); % xAM/Bmode, dB scale

% save data
clear pressure imMode;
if savedata
    save([saveName '_data_' datestr(now,'yymmdd-hh-MM-ss')],'P','saveName','sampROI_means','sampCNR','AM_Bmode_ratio','noiseROI_means','noiseROI_stds','voltage','PlateCoordinate','PlateSize','ROI_Centers','confScore','Nf');
end

%% plot ROI quants with microplateplot
% reshape ROI CNRs
% sampCNR new dimensions: well rows, well columns, frames, imaging modes
sampCNRs = permute(reshape(sampCNR, Nf, 2, PlateSize(2), PlateSize(1)), [4 3 1 2]);
AM_Bmode_ratios = permute(reshape(AM_Bmode_ratio, Nf, 1, PlateSize(2), PlateSize(1)), [4 3 1 2]);
confScore = permute(reshape(confScore, PlateSize(2), PlateSize(1)), [1 2]);

% find max signal achieved by each sample at any voltage
maxs_AM = squeeze(max(sampCNRs, [], 3));
maxs_AM_Bmode_ratio = squeeze(max(AM_Bmode_ratios, [], 3));

% make and save microplate plots
figure;
mpplot = microplateplot(maxs_AM(:,:,1));
colormap hot
colorbar
title('Max xAM signal achieved at any voltage')
mpplot;
savefig([saveName '_max-xAM'])

figure;
microplateplot(maxs_AM_Bmode_ratio(:,:))
colormap hot
colorbar
title('Max xAM:Bmode ratio signal achieved at any voltage')
mpplot;
savefig([saveName '_max-xAM-Bmode'])



%% function definitions
function updateNoiseROI(src,evt)
    disp(['ROI ' src.Tag ' moved. New position: ' mat2str(evt.CurrentPosition)])
    noiseMask = createMask(src); % Create mask for noise ROI
    
    % Calculate new quants
    wellIx = str2double(src.Tag(7:end));
    meanTemp = mean(nonzeros(evalin('base', sprintf('data(:,:,:, %d)', wellIx)).*noiseMask), [1 2 4], 'omitnan');
    stdTemp = std2(nonzeros(evalin('base', sprintf('data(:,:,:, %d)', wellIx)).*noiseMask));
    
    % Assign new quants into array
    evalin('base',sprintf('noiseROI_means(1,:, %d)=%d;',wellIx,meanTemp));
    disp(['New ' src.Tag ' mean: ' num2str(meanTemp)]);
    evalin('base',sprintf('noiseROI_stds(1,:, %d)=%d;',wellIx,stdTemp));
    disp(['New ' src.Tag ' STD: ' num2str(stdTemp)]);
end

function updateSampROI(src,evt)
    disp(['ROI ' src.Tag ' moved. New position: ' mat2str(evt.CurrentPosition)])
    sampMask = createMask(src); % Create mask for samp ROI
    
    % Calculate new quant
    wellIx = str2double(src.Tag(6:end));
    meanTemp = mean(nonzeros(evalin('base', sprintf('data(:,:, %d)', wellIx)).*sampMask), 'all', 'omitnan');
    
    % Assign new quant into array
    evalin('base',sprintf('sampROI_means(1, %d)=%d;',wellIx,meanTemp));
    disp(['New ' src.Tag ' mean: ' num2str(meanTemp)]);
end